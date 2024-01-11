'use strict';

export const RTCPeerConnection = (
    window.RTCPeerConnection ||
    window.webkitRTCPeerConnection ||
    window.mozRTCPeerConnection
).bind(window);

let peerConnection;
let streamId;
let sessionId;
let sessionClientAnswer;

let statsIntervalId;
let videoIsPlaying;
let lastBytesReceived;

const talkVideo = document.getElementById('talk-video');
talkVideo.setAttribute('playsinline', '');
const peerStatusLabel = document.getElementById('peer-status-label');
const iceStatusLabel = document.getElementById('ice-status-label');
const iceGatheringStatusLabel = document.getElementById('ice-gathering-status-label');
const signalingStatusLabel = document.getElementById('signaling-status-label');
const streamingStatusLabel = document.getElementById('streaming-status-label');



export async function talk(text) {
    // connectionState not supported in firefox
    var playing = sessionStorage.getItem('playing')
    if (playing == null) {
        playing = 'no'
    }

    if (text && playing == 'no') {
        sessionStorage.setItem('playing', 'yes')
        if (peerConnection?.signalingState === 'stable' || peerConnection?.iceConnectionState === 'connected') {
            const talkResponse = await fetchWithRetries(`api/createtalkstream`, {
                method: 'POST',
                headers: {
                    'Content-Type': 'application/json',
                },
                body: JSON.stringify({
                    text: text,
                    session_id: sessionId,
                    stream_id: streamId,
                }),
            });
            return talkResponse
        }
    } else if (playing == 'yes') {
        setTimeout(() => {
            talk(text)
        }, 1000)
    }

}

export async function connect() {
    if (peerConnection && peerConnection.connectionState === 'connected') {
        return;
    }

    stopAllStreams();
    closePC();

    const res = await fetch('api/createstream')

    const { id: newStreamId, offer, ice_servers: iceServers, session_id: newSessionId } = await res.json();

    streamId = newStreamId;
    sessionId = newSessionId;

    try {
        sessionClientAnswer = await createPeerConnection(offer, iceServers);
    } catch (e) {
        console.log('error during streaming setup', e);
        stopAllStreams();
        closePC();
        alert('An error occured during streaming setup');
        return;
    }

    const sdpResponse = await fetch(`api/startstream`, {
        method: 'POST',
        headers: {
            'Content-Type': 'application/json',
        },
        body: JSON.stringify({
            answer: sessionClientAnswer,
            session_id: sessionId,
            stream_id: streamId,
        }),
    });
    return
}

export async function destroyStream() {
    await fetch(`api/destroystream`, {
        method: 'POST',
        headers: {
            'Content-Type': 'application/json',
        },
        body: JSON.stringify({ session_id: sessionId, stream_id: streamId }),
    });

    stopAllStreams();
    closePC();
};

export function onIceGatheringStateChange() {
    iceGatheringStatusLabel.innerText = peerConnection.iceGatheringState;
    iceGatheringStatusLabel.className = 'iceGatheringState-' + peerConnection.iceGatheringState;
}
export function onIceCandidate(event) {
    console.log('onIceCandidate', event);
    if (event.candidate) {
        const { candidate, sdpMid, sdpMLineIndex } = event.candidate;

        fetch(`api/submitnetwork`, {
            method: 'POST',
            headers: {
                'Content-Type': 'application/json',
            },
            body: JSON.stringify({
                candidate,
                sdpMid,
                sdpMLineIndex,
                session_id: sessionId,
                stream_id: streamId,
            }),
        });
    }
}

export function onIceConnectionStateChange() {
    iceStatusLabel.innerText = peerConnection.iceConnectionState;
    iceStatusLabel.className = 'iceConnectionState-' + peerConnection.iceConnectionState;
    if (peerConnection.iceConnectionState === 'failed' || peerConnection.iceConnectionState === 'closed') {
        stopAllStreams();
        closePC();
    }
}
export function onConnectionStateChange() {
    // not supported in firefox
    peerStatusLabel.innerText = peerConnection.connectionState;
    peerStatusLabel.className = 'peerConnectionState-' + peerConnection.connectionState;
}
export function onSignalingStateChange() {
    signalingStatusLabel.innerText = peerConnection.signalingState;
    signalingStatusLabel.className = 'signalingState-' + peerConnection.signalingState;
}

export function onVideoStatusChange(videoIsPlaying, stream) {
    let status;
    if (videoIsPlaying) {
        sessionStorage.setItem('playing', 'yes')
        status = 'streaming';
        const remoteStream = stream;
        setVideoElement(remoteStream);
    } else {
        status = 'empty';
        playIdleVideo();
        sessionStorage.setItem('playing', 'no')
    }
    streamingStatusLabel.innerText = status;
    streamingStatusLabel.className = 'streamingState-' + status;
}

export function onTrack(event) {
    /**
     * The following code is designed to provide information about wether currently there is data
     * that's being streamed - It does so by periodically looking for changes in total stream data size
     *
     * This information in our case is used in order to show idle video while no talk is streaming.
     * To create this idle video use the POST https://api.d-id.com/talks endpoint with a silent audio file or a text script with only ssml breaks 
     * https://docs.aws.amazon.com/polly/latest/dg/supportedtags.html#break-tag
     * for seamless results use `config.fluent: true` and provide the same configuration as the streaming video
     */

    if (!event.track) return;

    statsIntervalId = setInterval(async () => {
        const stats = await peerConnection.getStats(event.track);
        stats.forEach((report) => {
            if (report.type === 'inbound-rtp' && report.mediaType === 'video') {
                const videoStatusChanged = videoIsPlaying !== report.bytesReceived > lastBytesReceived;

                if (videoStatusChanged) {
                    videoIsPlaying = report.bytesReceived > lastBytesReceived;
                    onVideoStatusChange(videoIsPlaying, event.streams[0]);
                }
                lastBytesReceived = report.bytesReceived;
            }
        });
    }, 500);
}

export async function createPeerConnection(offer, iceServers) {
    if (!peerConnection) {
        peerConnection = new RTCPeerConnection({ iceServers });
        peerConnection.addEventListener('icegatheringstatechange', onIceGatheringStateChange, true);
        peerConnection.addEventListener('icecandidate', onIceCandidate, true);
        peerConnection.addEventListener('iceconnectionstatechange', onIceConnectionStateChange, true);
        peerConnection.addEventListener('connectionstatechange', onConnectionStateChange, true);
        peerConnection.addEventListener('signalingstatechange', onSignalingStateChange, true);
        peerConnection.addEventListener('track', onTrack, true);
    }

    await peerConnection.setRemoteDescription(offer);
    console.log('set remote sdp OK');

    const sessionClientAnswer = await peerConnection.createAnswer();
    console.log('create local sdp OK');

    await peerConnection.setLocalDescription(sessionClientAnswer);
    console.log('set local sdp OK');

    return sessionClientAnswer;
}

export function setVideoElement(stream) {
    if (!stream) return;
    talkVideo.srcObject = stream;
    talkVideo.loop = false;

    // safari hotfix
    if (talkVideo.paused) {
        talkVideo
            .play()
            .then((_) => { })
            .catch((e) => { });
    }
}

export function playIdleVideo() {
    talkVideo.srcObject = undefined;
    //talkVideo.src = '../static/data/idleVideo.mp4';
    //talkVideo.loop = true;
}

export function stopAllStreams() {
    if (talkVideo.srcObject) {
        console.log('stopping video streams');
        talkVideo.srcObject.getTracks().forEach((track) => track.stop());
        talkVideo.srcObject = null;
    }
}

export function closePC(pc = peerConnection) {
    if (!pc) return;
    console.log('stopping peer connection');
    pc.close();
    pc.removeEventListener('icegatheringstatechange', onIceGatheringStateChange, true);
    pc.removeEventListener('icecandidate', onIceCandidate, true);
    pc.removeEventListener('iceconnectionstatechange', onIceConnectionStateChange, true);
    pc.removeEventListener('connectionstatechange', onConnectionStateChange, true);
    pc.removeEventListener('signalingstatechange', onSignalingStateChange, true);
    pc.removeEventListener('track', onTrack, true);
    clearInterval(statsIntervalId);
    iceGatheringStatusLabel.innerText = '';
    signalingStatusLabel.innerText = '';
    iceStatusLabel.innerText = '';
    peerStatusLabel.innerText = '';
    console.log('stopped peer connection');
    if (pc === peerConnection) {
        peerConnection = null;
    }
}

export const maxRetryCount = 3;
export const maxDelaySec = 4;

export async function fetchWithRetries(url, options, retries = 1) {
    try {
        return await fetch(url, options);
    } catch (err) {
        if (retries <= maxRetryCount) {
            const delay = Math.min(Math.pow(2, retries) / 4 + Math.random(), maxDelaySec) * 1000;

            await new Promise((resolve) => setTimeout(resolve, delay));

            console.log(`Request failed, retrying ${retries}/${maxRetryCount}. Error ${err}`);
            return fetchWithRetries(url, options, retries + 1);
        } else {
            throw new Error(`Max retries exceeded. error: ${err}`);
        }
    }
}