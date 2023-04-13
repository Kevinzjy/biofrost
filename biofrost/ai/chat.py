import json
import requests
from .utils import load_api_key


def chat(prompt, proxy_host='127.0.0.1', proxy_port='1087', api_key=None):
    """Chat with openAI GPT model

    Args:
        prompt (_type_): _description_
        proxy_host (str, optional): _description_. Defaults to '127.0.0.1'.
        proxy_port (str, optional): _description_. Defaults to '1087'.
        api_key (_type_, optional): _description_. Defaults to None.

    Returns:
        _type_: _description_
    """
    proxies = {
        "http": f"http://{proxy_host}:{proxy_port}",
        "https": f"http://{proxy_host}:{proxy_port}",
    }

    if api_key is None:
        api_key = load_api_key()
    headers = {
        'Content-Type': 'application/json',
        'Authorization': f'Bearer {api_key}',
    }

    url = 'https://api.openai.com/v1/chat/completions'
    model_name = 'gpt-3.5-turbo'
    json_data = {'model': model_name, 'messages': [{'role': 'user', 'content': prompt}]}

    response = requests.post(url, headers=headers, json=json_data, proxies=proxies)
    ret = json.loads(response.text.strip())
    text = ret['choices'][0]['message']['content']
    return text