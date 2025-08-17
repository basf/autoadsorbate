"""Site utility functions."""

import pandas as pd


def make_site_info_writable(site_info):
    """
    Convert site information to a writable format.
    
    Args:
        site_info: Site information dictionary or object
        
    Returns:
        dict: Writable site information
    """
    if isinstance(site_info, dict):
        return site_info.copy()
    elif hasattr(site_info, '__dict__'):
        return site_info.__dict__.copy()
    else:
        return {'site_info': str(site_info)}
