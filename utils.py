import re     # Generate cartesian products (w-mers)
import random # Generate random sequence q


def find(substring, string):
    '''
    Return starting indices of all occurences
    of substring in string with overlap
    '''
    return tuple(
        match.start()
        for match in re.finditer(f'(?={substring})', string)
    )


def random_str(charset, length):
    '''
    Generate random string from charset of fixed length
    '''
    return ''.join(random.choice(charset) for _ in range(length))