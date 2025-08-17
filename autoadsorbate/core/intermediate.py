"""Intermediate class for handling reaction intermediates."""


class Intermediate:
    """
    Base class for initializing reaction intermediates.

    Attributes:
        ActiveSite: The active site for the intermediate.
        fragments: A list of fragments associated with the intermediate.
    """

    def __init__(self, ActiveSite, fragments=None):
        """
        Initialize attributes.

        Args:
            ActiveSite: The active site for the intermediate.
            fragments (list, optional): A list of fragments associated with the intermediate. Defaults to an empty list.
        """
        self.ActiveSite = ActiveSite
        self.fragments = fragments if fragments is not None else []
