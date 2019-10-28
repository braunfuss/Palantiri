Common errors and fixes
==============

A list of common errors and how to avoid/fix them.

1. Time window too short
---------------
If the time window given by duration and forerun is too short in comparison to the
step and window size, no semblance can be reliably calculated.
The following error will appear:


.. highlight:: console

::

    $ ValueError: zero-size array to reduction operation maximum which has no identity

A simple fix is to choose at least an forerun+duration length equal to the step size.
