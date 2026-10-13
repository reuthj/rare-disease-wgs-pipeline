"""Interpret the main variant classification without resolving source conflicts."""
import re


def classification_tokens(value):
    text = str(value).strip().casefold().replace('_', ' ')
    return [re.sub(r'\s+', ' ', part.strip()) for part in
            re.split(r'[|/;,+]', text) if part.strip()]


def is_main_plp(value):
    tokens = classification_tokens(value)
    return bool(tokens) and all(token in {'p', 'lp', 'pathogenic', 'likely pathogenic'}
                                for token in tokens)
