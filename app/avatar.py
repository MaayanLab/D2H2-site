import requests
import os
from dotenv import load_dotenv
load_dotenv()

from functools import lru_cache


D_ID_API_KEY = os.environ.get('D_ID_API_KEY')
SOURCE_URL = os.environ.get('SOURCE_URL')
def create_stream():
    """Starts the stream of the avatar"""
    url = "https://api.d-id.com/talks/streams"
    

    payload = {
        "source_url": SOURCE_URL,
    }

    headers = {
        "accept": "application/json",
        "content-type": "application/json",
        "Authorization": f"Basic {D_ID_API_KEY}",
    }

    response = requests.post(url, json=payload, headers=headers)
    return response.json()

def start_stream(stream_id: str, session_id: str, answer: dict):
    url = f"https://api.d-id.com/talks/streams/{stream_id}/sdp"

    payload = {
        "session_id": session_id,
        "answer": answer
    }
    headers = {
        "accept": "application/json",
        "content-type": "application/json",
        "Authorization": f"Basic {D_ID_API_KEY}"
    }

    response = requests.post(url, json=payload, headers=headers)
    return response.json()


def submit_network(stream_id: str, session_id: str, candidate, sdpMid, sdpMLineIndex):

    url = f"https://api.d-id.com/talks/streams/{stream_id}/ice"

    payload = {
        "candidate": candidate,
        "sdpMid": sdpMid,
        "sdpMLineIndex": sdpMLineIndex,
        "session_id": session_id
    }
    headers = {
        "accept": "application/json",
        "content-type": "application/json",
        "Authorization": f"Basic {D_ID_API_KEY}"
    }

    response = requests.post(url, json=payload, headers=headers)

    return response.json()

@lru_cache()
def create_talk_stream(stream_id, session_id, text):

    url = f"https://api.d-id.com/talks/streams/{stream_id}"
    print(stream_id, session_id, text)
    payload = {
        "script": {
            "type": "text",
            "subtitles": "false",
            "provider": {
                "type": "microsoft",
                "voice_id": "en-US-AndrewNeural"
            },
            "ssml": False,
            "input": text 
        },
        "config": {
            "fluent": True,
            "pad_audio": "0.0",
            "align_driver": True,
            "auto_match": True
        },
        "session_id": session_id,
    }
    headers = {
        "content-type": "application/json",
        "Authorization": f"Basic {D_ID_API_KEY}"
    }

    response = requests.post(url, json=payload, headers=headers)
    return response.json()


def destroy_stream(stream_id, seesion_id):

    url = f"https://api.d-id.com/talks/streams/{stream_id}"

    payload = { "session_id": seesion_id }
    headers = {
        "accept": "application/json",
        "content-type": "application/json",
        "Authorization": f"Basic {D_ID_API_KEY}"
    }

    response = requests.delete(url, json=payload, headers=headers)
    return response


def generate_idle_video():

    url = "https://api.d-id.com/talks"

    payload = {
        "script": {
            "type": "text",
            "subtitles": False,
            "provider": {
                "type": "microsoft",
                "voice_id": "en-US-JennyNeural",
            },
            "ssml": True,
            "input": "<break time=\"5s\"/>"
        },
        "config": {
            "fluent": True,
            "pad_audio": 0,
            "align_driver": True,
            "auto_match": True,
        },
        "source_url": SOURCE_URL,
        "webhook": "https://host.domain.tld/to/webhook",
        "persist": True
    }
    headers = {
        "accept": "application/json",
        "content-type": "application/json",
        "Authorization": f"Basic {D_ID_API_KEY}"
    }

    response = requests.post(url, json=payload, headers=headers)