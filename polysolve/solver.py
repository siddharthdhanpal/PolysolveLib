from cmath import sqrt
import math

from cowsay import cow

CBRT_UNITY_IM = sqrt(3)/2 * 1j

def quadratic(a: float, b: float, c: float) -> tuple[float, float]:
    """
    Solves the roots of a quadratic equation.

    Uses the quadratic formula. Result must be real.

    Parameters
    ----------
    a
       :math:`x^2` coefficient.
    b
       :math:`x` coefficient.
    c
       Constant value.

    Returns
    -------
    tuple[float, float]
        Positive and negative roots of quadratic.

    Raises
    ------
    ValueError
        Discriminant < 0 implying imaginary root.

    Notes
    -----
    Equation of the form:

    .. math::

        ax^{2} + bx + c

    Examples
    --------
    >>> quadratic(1, 2, 0)
    (0.0, -2.0)
    >>> quadratic(3., 0., -1.)
    (0.5773502691896257, -0.5773502691896257)

    See Also
    --------
    numpy.polyval : Evaluate polynomial at point.

    References
    ----------
    .. [1] O. McNoleg, "The integration of GIS, remote sensing,
           expert systems ...
    """

    det = b**2 - (4*a*c)

    return ((-b + sqrt(det)) / (2*a),
            (-b - sqrt(det)) / (2*a))


# def quadratic(a, b, c):
#     det = b**2 - (4*a*c)

#     if math.isclose(det, 0):
#         cow("Degenerate MOOoo-ts")

#     return ((-b + sqrt(det)) / (2*a), (-b - sqrt(det)) / (2*a))

def cubic(a: float, b: float, c: float, d: float) -> tuple[float, float, float]:
    """
    Solves the roots of a quadratic equation.

    Uses the quadratic formula. Result must be real.

    Parameters
    ----------
    a
       :math:`x^3` coefficient.
    b
       :math:`x^2` coefficient.
    c
       :math:`x` coefficient.    
    d
       :Constant value.

    Returns
    -------
    tuple[float, float, float]
        Positive and negative roots of cubic.

    Raises
    ------
    ValueError
        Discriminant < 0 implying imaginary root.

    Notes
    -----
    Equation of the form:

    .. math::

        ax^{3} + bx^{2} + cx+d

    Examples
    --------
    >>> quadratic(1, 2, 0)
    (0.0, -2.0)
    >>> quadratic(3., 0., -1.)
    (0.5773502691896257, -0.5773502691896257)

    See Also
    --------
    numpy.polyval : Evaluate polynomial at point.

    References
    ----------
    .. [1] O. McNoleg, "The integration of GIS, remote sensing,
           expert systems ...
    """

    q = (3*a*c - b**2) / (9*a**2)
    r = (9*a*b*c - 27*a**2*d - 2*b**3) / (54*a**3)

    s = (r + sqrt(q**3 + r**2))**(1/3)
    t = (r - sqrt(q**3 + r**2))**(1/3)

    x1 = s + t - (b/3*a)
    x2 = -(s + t)/2 - (b/3*a) + CBRT_UNITY_IM * (s - t)
    x3 = -(s + t)/2 - (b/3*a) - CBRT_UNITY_IM * (s - t)

    if any(x == x1 for x in (x2, x3)):
        cow("Degenerate MOOoo-ts")

    return (x1, x2, x3)
