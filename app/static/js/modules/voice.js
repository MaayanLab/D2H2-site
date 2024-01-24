
const recorder = new MicRecorder({
    bitRate: 128
});

const ctx = new AudioContext();


export async function play_message(message) {
    var playing = sessionStorage.getItem('playing')
    if (playing == null) {
        playing = 'no'
    }

    if (message && playing == 'no') {
        sessionStorage.setItem('playing', 'yes')
        const response = await fetch("/api/speak_message", {
            method: 'POST',
            body: JSON.stringify({
                'text':
                    message
            }),
            headers: {
                'Content-Type': 'application/json',
            }
        })
        const playSound = ctx.createBufferSource();
        const arrayBuffer = await response.arrayBuffer()
        const audio = await ctx.decodeAudioData(arrayBuffer)

        playSound.buffer = audio;
        playSound.connect(ctx.destination);
        playSound.start(ctx.currentTime);
        setTimeout(() => {
            sessionStorage.setItem('playing', 'no')
        }, audio.duration * 1000)
    } else if (playing == 'yes') {
        setTimeout(() => {
            play_message(message)
        }, 1000)
    }
}

    export async function start_recording() {
        try {
            console.log("recording")
            recorder.start()
        } catch (error) {
            console.error(error)
        }
    }


    export async function stop_recording() {
        try {
            console.log("saving")
            const [buffer, blob] = await recorder.stop().getMp3()
            const formData = new FormData();
            formData.append("file", blob);
            console.log(blob)
            if (blob.size < 10) {
                return {}
            }
            const response = await fetch("api/transcribe_message",
                {
                    method: 'POST',
                    body: formData,
                    headers: {
                        'Accept': 'application/json',
                    }
                });
            const transcribed = await response.json()
            console.log(transcribed)
            return transcribed
        } catch (error) {
            console.error(error)
        }
    }