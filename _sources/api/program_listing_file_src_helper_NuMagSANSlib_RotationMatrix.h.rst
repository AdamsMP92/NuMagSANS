
.. _program_listing_file_src_helper_NuMagSANSlib_RotationMatrix.h:

Program Listing for File NuMagSANSlib_RotationMatrix.h
======================================================

|exhale_lsh| :ref:`Return to documentation for file <file_src_helper_NuMagSANSlib_RotationMatrix.h>` (``src/helper/NuMagSANSlib_RotationMatrix.h``)

.. |exhale_lsh| unicode:: U+021B0 .. UPWARDS ARROW WITH TIP LEFTWARDS

.. code-block:: cpp

   #pragma once
   
   #include <algorithm>
   #include <array>
   #include <cctype>
   #include <cmath>
   #include <stdexcept>
   #include <string>
   
   inline std::string NormalizeEulerConvention(std::string convention) {
       std::transform(convention.begin(), convention.end(), convention.begin(),
                      [](unsigned char character) { return static_cast<char>(std::tolower(character)); });
       return convention;
   }
   
   inline bool IsSupportedEulerConvention(const std::string& convention) {
       static constexpr std::array<const char*, 12> supported = {
           "xyx", "xzx", "yxy", "yzy", "zxz", "zyz", "xyz", "xzy", "yxz", "yzx", "zxy", "zyx",
       };
   
       const std::string normalized = NormalizeEulerConvention(convention);
       return std::find_if(supported.begin(), supported.end(),
                           [&](const char* candidate) { return normalized == candidate; }) != supported.end();
   }
   
   inline void RotationMatrix_x(float angle, float* rotation_matrix) {
       const float cosine = std::cos(angle);
       const float sine = std::sin(angle);
   
       // Column-major storage: R(row, column) = R[row + 3 * column].
       rotation_matrix[0] = 1.0f;
       rotation_matrix[1] = 0.0f;
       rotation_matrix[2] = 0.0f;
   
       rotation_matrix[3] = 0.0f;
       rotation_matrix[4] = cosine;
       rotation_matrix[5] = sine;
   
       rotation_matrix[6] = 0.0f;
       rotation_matrix[7] = -sine;
       rotation_matrix[8] = cosine;
   }
   
   inline void RotationMatrix_y(float angle, float* rotation_matrix) {
       const float cosine = std::cos(angle);
       const float sine = std::sin(angle);
   
       rotation_matrix[0] = cosine;
       rotation_matrix[1] = 0.0f;
       rotation_matrix[2] = -sine;
   
       rotation_matrix[3] = 0.0f;
       rotation_matrix[4] = 1.0f;
       rotation_matrix[5] = 0.0f;
   
       rotation_matrix[6] = sine;
       rotation_matrix[7] = 0.0f;
       rotation_matrix[8] = cosine;
   }
   
   inline void RotationMatrix_z(float angle, float* rotation_matrix) {
       const float cosine = std::cos(angle);
       const float sine = std::sin(angle);
   
       rotation_matrix[0] = cosine;
       rotation_matrix[1] = sine;
       rotation_matrix[2] = 0.0f;
   
       rotation_matrix[3] = -sine;
       rotation_matrix[4] = cosine;
       rotation_matrix[5] = 0.0f;
   
       rotation_matrix[6] = 0.0f;
       rotation_matrix[7] = 0.0f;
       rotation_matrix[8] = 1.0f;
   }
   
   inline void RotationMatrix_axis(char axis, float angle, float* rotation_matrix) {
       switch (static_cast<char>(std::tolower(static_cast<unsigned char>(axis)))) {
       case 'x':
           RotationMatrix_x(angle, rotation_matrix);
           return;
       case 'y':
           RotationMatrix_y(angle, rotation_matrix);
           return;
       case 'z':
           RotationMatrix_z(angle, rotation_matrix);
           return;
       default:
           throw std::invalid_argument("Euler rotation axis must be x, y, or z.");
       }
   }
   
   inline void LeftMultiply3x3(float* left, const float* multiplier) {
       // Updates left in place as left <- multiplier * left. Matrices use
       // column-major storage: R(row, column) = R[row + 3 * column].
       float product[9];
   
       for (int column = 0; column < 3; ++column) {
           for (int row = 0; row < 3; ++row) {
               product[row + 3 * column] = multiplier[row + 3 * 0] * left[0 + 3 * column] +
                                           multiplier[row + 3 * 1] * left[1 + 3 * column] +
                                           multiplier[row + 3 * 2] * left[2 + 3 * column];
           }
       }
   
       for (int index = 0; index < 9; ++index) {
           left[index] = product[index];
       }
   }
   
   inline void Multiply_RotmatEuler_3x3(float angle_1, float angle_2, float angle_3, const std::string& convention,
                                        float* rotation_matrix) {
       const std::string normalized = NormalizeEulerConvention(convention);
       if (!IsSupportedEulerConvention(normalized)) {
           throw std::invalid_argument("Unsupported Euler rotation convention: " + convention);
       }
   
       float first_axis_rotation[9];
       float second_axis_rotation[9];
       float third_axis_rotation[9];
   
       RotationMatrix_axis(normalized[0], angle_1, first_axis_rotation);
       RotationMatrix_axis(normalized[1], angle_2, second_axis_rotation);
       RotationMatrix_axis(normalized[2], angle_3, third_axis_rotation);
   
       // Starting from identity, this produces
       // R = R_axis1(angle_1) * R_axis2(angle_2) * R_axis3(angle_3).
       LeftMultiply3x3(rotation_matrix, third_axis_rotation);
       LeftMultiply3x3(rotation_matrix, second_axis_rotation);
       LeftMultiply3x3(rotation_matrix, first_axis_rotation);
   }
   
   inline void ComputeEulerRotationMatrixDegrees_3x3(float angle_1_degrees, float angle_2_degrees, float angle_3_degrees,
                                                     const std::string& convention, float* rotation_matrix) {
       constexpr float degrees_to_radians = 3.14159265358979323846f / 180.0f;
   
       // Global RotMat angles use degrees, whereas object-wise RotData angles use
       // radians. Both interfaces share the same Euler convention and angle order.
       rotation_matrix[0] = 1.0f;
       rotation_matrix[1] = 0.0f;
       rotation_matrix[2] = 0.0f;
       rotation_matrix[3] = 0.0f;
       rotation_matrix[4] = 1.0f;
       rotation_matrix[5] = 0.0f;
       rotation_matrix[6] = 0.0f;
       rotation_matrix[7] = 0.0f;
       rotation_matrix[8] = 1.0f;
   
       Multiply_RotmatEuler_3x3(angle_1_degrees * degrees_to_radians, angle_2_degrees * degrees_to_radians,
                                angle_3_degrees * degrees_to_radians, convention, rotation_matrix);
   }
