"""
Module for working with cross-section properties of beam elements.
"""
module CrossSectionModule

using Dierckx
using LinearAlgebra: norm
using FinEtools

abstract type AbstractCrossSectionType end

struct CrossSectionCircle{F} <: AbstractCrossSectionType
    shape::String
    # Function to return A, J, I1, I2, I3
    #     `A` = cross-sectional area of the rod
    #     `I1`= moment of inertia of the cross-section, twisting about the beam axis 
    #     `I2`, `I3`=moments of inertia, bending about the x2-axis (I2) or x3-axis (I3)
    #     `J`= torsion constant, (to give torsion stiffness =G*J)
    #     `A2s` = effective area for shear in the direction of x2
    #     `A3s` = effective area for shear in the direction of x3
    #     x1x2_vector = vector to determine the local x1--x2 plane
    # The parameter is the normalized arc length.
    # A, I1, I2, I3, J, A2s, A3s, x1x2_vector = parameters(s)
    parameters::F
    label::Int
end

"""
    CrossSectionCircle(radius, x1x2_vector; label = 0)

Create a circular section for use with a Bernoulli model.
"""
function CrossSectionCircle(radius, x1x2_vector; label = 0)
    function parameters(s)
        R = radius(s)
        A = pi * R^2
        I1 = pi / 2 * R^4
        I2 = pi / 4 * R^4
        I3 = pi / 4 * R^4
        J = pi / 2 * R^4
        A2s = Inf # Bernoulli is the default
        A3s = Inf # Bernoulli is the default
        return (
            A = A,
            I1 = I1,
            I2 = I2,
            I3 = I3,
            J = J,
            A2s = A2s,
            A3s = A3s,
            x1x2_vector = x1x2_vector(s),
            dimensions = [R],
        )
    end
    return CrossSectionCircle("circle", parameters, label)
end

"""
    CrossSectionCircle(radius, x1x2_vector, shear_resisting_area_factor = 6/7; label = 0)

Create a circular section for use with a Timoshenko model.

Commonly accepted value is `shear_resisting_area_factor = 6/7`.
"""
function CrossSectionCircle(radius, x1x2_vector, shear_resisting_area_factor; label = 0)
    function parameters(s)
        R = radius(s)
        A = pi * R^2
        I1 = pi / 2 * R^4
        I2 = pi / 4 * R^4
        I3 = pi / 4 * R^4
        J = pi / 2 * R^4
        A2s = A * shear_resisting_area_factor # Timoshenko
        A3s = A * shear_resisting_area_factor # Timoshenko
        return (
            A = A,
            I1 = I1,
            I2 = I2,
            I3 = I3,
            J = J,
            A2s = A2s,
            A3s = A3s,
            x1x2_vector = x1x2_vector(s),
            dimensions = [R],
        )
    end
    return CrossSectionCircle("circle", parameters, label)
end

struct CrossSectionHollowCircle{F} <: AbstractCrossSectionType
    shape::String
    parameters::F
    label::Int
end

function CrossSectionHollowCircle(innerradius, outerradius, x1x2_vector; label = 0)
    function parameters(s)
        Rext = outerradius(s)
        Rint = innerradius(s)
        A = pi * (Rext^2 - Rint^2)
        I2 = pi / 4 * (Rext^4 - Rint^4)
        I3 = pi / 4 * (Rext^4 - Rint^4)
        I1 = I2 + I3
        J = pi / 2 * (Rext^4 - Rint^4)
        A2s = Inf # Bernoulli is the default
        A3s = Inf # Bernoulli is the default
        return (
            A = A,
            I1 = I1,
            I2 = I2,
            I3 = I3,
            J = J,
            A2s = A2s,
            A3s = A3s,
            x1x2_vector = x1x2_vector(s),
            dimensions = [Rext, Rint],
        )
    end
    return CrossSectionHollowCircle("hollow circle", parameters, label)
end

struct CrossSectionRectangle{F} <: AbstractCrossSectionType
    shape::String
    parameters::F
    label::Int
end

"""
    CrossSectionRectangle(d2f, d3f, x1x2_vector; label = 0)

Create a rectangular section for use with a Bernoulli model.
"""
function CrossSectionRectangle(d2f, d3f, x1x2_vector; label = 0)
    function parameters(s)
        d2 = d2f(s)
        d3 = d3f(s)
        a = max(d2, d3)
        b = min(d2, d3)
        itp = Spline1D(
            [1, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 6.0, 10, 20, 40, 80, 200, 2000],
            [
                0.141,
                0.196,
                0.229,
                0.249,
                0.263,
                0.281,
                0.291,
                0.299,
                0.312,
                0.317,
                0.325,
                0.33,
                1 / 3,
                1 / 3,
            ],
        )
        c = itp(a / b)
        A = d2 * d3
        I2 = d2 * d3^3 / 12
        I3 = d2^3 * d3 / 12
        I1 = I2 + I3
        J = c * a * b^3
        A2s = Inf # Bernoulli is the default
        A3s = Inf # Bernoulli is the default
        return (
            A = A,
            I1 = I1,
            I2 = I2,
            I3 = I3,
            J = J,
            A2s = A2s,
            A3s = A3s,
            x1x2_vector = x1x2_vector(s),
            dimensions = [d2, d3],
        )
    end
    return CrossSectionRectangle("rectangle", parameters, label)
end

"""
    CrossSectionRectangle(d2f, d3f, x1x2_vector; label = 0)

Create a rectangular section for use with a Timoshenko model.

Commonly accepted value is `shear_resisting_area_factor = 5/6`.
"""
function CrossSectionRectangle(
    d2f,
    d3f,
    x1x2_vector,
    shear_resisting_area_factor;
    label = 0,
)
    function parameters(s)
        d2 = d2f(s)
        d3 = d3f(s)
        a = max(d2, d3)
        b = min(d2, d3)
        itp = Spline1D(
            [1, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 6.0, 10, 20, 40, 80, 200, 2000],
            [
                0.141,
                0.196,
                0.229,
                0.249,
                0.263,
                0.281,
                0.291,
                0.299,
                0.312,
                0.317,
                0.325,
                0.33,
                1 / 3,
                1 / 3,
            ],
        )
        c = itp(a / b)
        A = d2 * d3
        I2 = d2 * d3^3 / 12
        I3 = d2^3 * d3 / 12
        I1 = I2 + I3
        J = c * a * b^3
        A2s = A * shear_resisting_area_factor # Timoshenko 
        A3s = A * shear_resisting_area_factor # Timoshenko
        return (
            A = A,
            I1 = I1,
            I2 = I2,
            I3 = I3,
            J = J,
            A2s = A2s,
            A3s = A3s,
            x1x2_vector = x1x2_vector(s),
            dimensions = [d2, d3],
        )
    end
    return CrossSectionRectangle("rectangle", parameters, label)
end

end # module
