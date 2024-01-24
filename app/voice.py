import os
from pathlib import Path
from openai import OpenAI
import io

exports = dict(
  openai=OpenAI(),
) if 'OPENAI_API_KEY' in os.environ else dict()


async def speak(text):
    if ('openai' in exports):
        client = exports["openai"]
        response = client.audio.speech.create(
            model="tts-1",
            voice="alloy",
            input=text
        )
        return response.content
    else:
        raise {'status_code': 400, 'detail': "Invalid Input"}

async def transcribe(file):
    if ('openai' in exports):
        audio_file = io.BytesIO(file.read())
        audio_file.name = "speech.mp3"
        client = exports["openai"]
        transcript = client.audio.transcriptions.create(
            model="whisper-1", 
            file=audio_file
        )
        return transcript.text
    else:
        raise {'status_code': 400, 'detail': "Invalid Input"}
