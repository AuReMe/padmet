import os


def is_valid_path(filepath):
    """Return True if filepath is valid
    
    Args:
        filepath (str): path to file
    
    Returns:
        bool: True if path exists, False otherwise
    """
    if filepath and not os.access(filepath, os.W_OK):
        try:
            open(filepath, 'w').close()
            os.unlink(filepath)
            return True
        except OSError:
            return False
    else:  # path is accessible
        return True
