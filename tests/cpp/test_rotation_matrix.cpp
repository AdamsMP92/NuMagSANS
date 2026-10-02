#include "helper/NuMagSANSlib_RotationMatrix.h"

#include <array>
#include <cmath>
#include <stdexcept>
#include <string>

namespace {

constexpr float pi = 3.14159265358979323846f;
constexpr float tolerance = 1.0e-6f;

void require(bool condition, const std::string& message) {
    if (!condition) {
        throw std::runtime_error(message);
    }
}

std::array<float, 9> identity_matrix() {
    return {1.0f, 0.0f, 0.0f, 0.0f, 1.0f, 0.0f, 0.0f, 0.0f, 1.0f};
}

std::array<float, 9> legacy_zyz_matrix(float alpha, float beta, float gamma) {
    const float ca = std::cos(alpha);
    const float sa = std::sin(alpha);
    const float cb = std::cos(beta);
    const float sb = std::sin(beta);
    const float cg = std::cos(gamma);
    const float sg = std::sin(gamma);

    // Closed-form equivalent of the historical
    // Rz(alpha) * Ry(beta) * Rz(gamma) implementation, stored column-major.
    return {
        ca * cb * cg - sa * sg,
        sa * cb * cg + ca * sg,
        -sb * cg,
        -ca * cb * sg - sa * cg,
        -sa * cb * sg + ca * cg,
        sb * sg,
        ca * sb,
        sa * sb,
        cb,
    };
}

void assert_matrix_close(const std::array<float, 9>& actual, const std::array<float, 9>& expected) {
    for (std::size_t index = 0; index < actual.size(); ++index) {
        require(std::abs(actual[index] - expected[index]) < tolerance,
                "Rotation matrices differ at column-major index " + std::to_string(index));
    }
}

void test_supported_conventions() {
    const std::array<std::string, 12> supported = {
        "xyx", "xzx", "yxy", "yzy", "zxz", "zyz", "xyz", "xzy", "yxz", "yzx", "zxy", "zyx",
    };

    for (const std::string& convention : supported) {
        require(IsSupportedEulerConvention(convention), "Expected supported convention: " + convention);
        std::array<float, 9> rotation = identity_matrix();
        Multiply_RotmatEuler_3x3(0.0f, 0.0f, 0.0f, convention, rotation.data());
        assert_matrix_close(rotation, identity_matrix());
    }

    require(IsSupportedEulerConvention("XYZ"), "Convention matching should be case-insensitive");
    require(!IsSupportedEulerConvention("xxy"), "xxy must not be accepted");
    require(!IsSupportedEulerConvention("abc"), "abc must not be accepted");
}

void test_xyz_x_rotation() {
    std::array<float, 9> rotation = identity_matrix();
    Multiply_RotmatEuler_3x3(pi / 2.0f, 0.0f, 0.0f, "xyz", rotation.data());

    const std::array<float, 9> expected = {
        1.0f, 0.0f, 0.0f, 0.0f, 0.0f, 1.0f, 0.0f, -1.0f, 0.0f,
    };
    assert_matrix_close(rotation, expected);
}

void test_zyz_default_matches_equivalent_x_rotation() {
    std::array<float, 9> xyz_rotation = identity_matrix();
    std::array<float, 9> zyz_rotation = identity_matrix();

    Multiply_RotmatEuler_3x3(pi / 2.0f, 0.0f, 0.0f, "xyz", xyz_rotation.data());
    Multiply_RotmatEuler_3x3(-pi / 2.0f, pi / 2.0f, pi / 2.0f, "zyz", zyz_rotation.data());

    assert_matrix_close(zyz_rotation, xyz_rotation);
}

void test_zyz_matches_historical_formula() {
    constexpr float alpha = 0.31f;
    constexpr float beta = -0.47f;
    constexpr float gamma = 1.02f;
    std::array<float, 9> rotation = identity_matrix();

    Multiply_RotmatEuler_3x3(alpha, beta, gamma, "zyz", rotation.data());

    assert_matrix_close(rotation, legacy_zyz_matrix(alpha, beta, gamma));
}

void test_global_degree_angles_follow_same_zyz_order() {
    constexpr float alpha_degrees = 17.0f;
    constexpr float beta_degrees = -31.0f;
    constexpr float gamma_degrees = 49.0f;
    constexpr float degrees_to_radians = pi / 180.0f;

    std::array<float, 9> global_rotation{};
    ComputeEulerRotationMatrixDegrees_3x3(alpha_degrees, beta_degrees, gamma_degrees, "zyz", global_rotation.data());

    const std::array<float, 9> expected = legacy_zyz_matrix(
        alpha_degrees * degrees_to_radians, beta_degrees * degrees_to_radians, gamma_degrees * degrees_to_radians);
    assert_matrix_close(global_rotation, expected);
}

void test_global_xyz_first_angle_is_x_rotation() {
    std::array<float, 9> rotation{};
    ComputeEulerRotationMatrixDegrees_3x3(90.0f, 0.0f, 0.0f, "xyz", rotation.data());

    const std::array<float, 9> expected = {
        1.0f, 0.0f, 0.0f, 0.0f, 0.0f, 1.0f, 0.0f, -1.0f, 0.0f,
    };
    assert_matrix_close(rotation, expected);
}

void test_convention_controls_multiplication_order() {
    std::array<float, 9> rotation = identity_matrix();
    Multiply_RotmatEuler_3x3(pi / 2.0f, pi / 2.0f, 0.0f, "xyz", rotation.data());

    // R = Rx(pi/2) * Ry(pi/2), stored column-major.
    const std::array<float, 9> expected = {
        0.0f, 1.0f, 0.0f, 0.0f, 0.0f, 1.0f, 1.0f, 0.0f, 0.0f,
    };
    assert_matrix_close(rotation, expected);
}

void test_invalid_convention_is_rejected() {
    std::array<float, 9> rotation = identity_matrix();
    bool rejected = false;

    try {
        Multiply_RotmatEuler_3x3(0.0f, 0.0f, 0.0f, "xxy", rotation.data());
    } catch (const std::invalid_argument&) {
        rejected = true;
    }

    require(rejected, "Invalid convention must raise std::invalid_argument");
}

} // namespace

int main() {
    test_supported_conventions();
    test_xyz_x_rotation();
    test_zyz_default_matches_equivalent_x_rotation();
    test_zyz_matches_historical_formula();
    test_global_degree_angles_follow_same_zyz_order();
    test_global_xyz_first_angle_is_x_rotation();
    test_convention_controls_multiplication_order();
    test_invalid_convention_is_rejected();
    return 0;
}
