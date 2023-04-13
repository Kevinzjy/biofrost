from pathlib import Path


def load_api_key():
    key_file = Path.home() / ".openai_key"
    if key_file.exists():
        with open(key_file, 'r') as f:
            return f.readline().rstrip()
    else:
        return KeyError("No openai API key found. please save it in ~/.openai_key.pub.")