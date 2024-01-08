
const recorder = new MicRecorder({
    bitRate: 128
});

export async function play_message(message) {
    if (message) {
        const response = await fetch("/api/speak_message", {
            method: 'POST',
            body: JSON.stringify({'text':
                message
            }),
            headers: {
                'Content-Type': 'application/json',
            }
        })
        const ctx = new AudioContext();
        const arrayBuffer = await response.arrayBuffer()
        const audio = await ctx.decodeAudioData(arrayBuffer)
        const playSound = ctx.createBufferSource();
        playSound.buffer = audio;
        playSound.connect(ctx.destination);
        playSound.start(ctx.currentTime);
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