#include "triangulation.h"
#include "matrix_algo.h"
#include <easy3d/optimizer/optimizer_lm.h>

using namespace easy3d;

/** Minimise using SVD (from matrix_algo.h Easy3D) and return the last column of V
 * Description: This function takes a matrix M and performs SVD decomposition on it.
 * input: Matrix M
 * output: Vector - last column of V
 */
Vector minimize_using_svd(Matrix M) {
    int num_rows = M.rows();
    int num_cols = M.cols();
    Matrix U(num_rows, num_rows, 0.0);
    Matrix D(num_rows, num_cols, 0.0);
    Matrix V(num_cols, num_cols, 0.0);

    svd_decompose(M, U, D, V);

    // Print the matrix
//    std::cout << "Matrix v: " << V.get_column(V.cols() - 1) << std::endl;
//    std::cout << "-----------------" << std::endl;

    return V.get_column(V.cols() - 1);
}

/** Create M
 * Description: This function creates the M matrix
 * Input: Matrix K_i - calibration matrix, Matrix R_i - rotation matrix, Vector3D t_i - translation vector
 * Output: Matrix M_i - projection matrix
 */
Matrix createM(Matrix K_i, Matrix R_i, Vector3D t_i) {

    // Inistialise M matrix 3x4 with 0.0
    Matrix M_i(3, 4, 0.0);

    M_i.set_column(0, R_i.get_column(0));
    M_i.set_column(1, R_i.get_column(1));
    M_i.set_column(2, R_i.get_column(2));

    M_i.set_column(3, t_i);

    // Apply the calibration matrix
    M_i = K_i * M_i;

    // Uncomment to print the matrix
//    std::cout << "Matrix M: " << M_i << std::endl;

    return M_i;
}


/** Create A
 * Description: This function creates the A matrix
 * Input: Vector2D p0 - 2D image point in the first image, Vector2D p1 - 2D image point in the second image,
 *          Matrix M_i - projection matrix of the second camera, Matrix M_l - projection matrix of the first camera
 * Output: Matrix A - 4x4 matrix
 */
Matrix createA(Vector2D p0, Vector2D p1, Matrix M_i, Matrix M_l) {

    // Inistialise A matrix 4x4 with 0.0
    Matrix A(4, 4, 0.0);

    A.set_row(0, p0[0] * M_l.transpose().get_column(2) - M_l.transpose().get_column(0));
    A.set_row(1, p0[1] * M_l.transpose().get_column(2) - M_l.transpose().get_column(1));
    A.set_row(2, p1[0] * M_i.transpose().get_column(2) - M_i.transpose().get_column(0));
    A.set_row(3, p1[1] * M_i.transpose().get_column(2) - M_i.transpose().get_column(1));

    return A;
}

/**
 * @return True on success, otherwise false. On success, the reconstructed 3D points must be written to 'points_3d'
 *      and the recovered relative pose must be written to R and t.
 */

bool Triangulation::triangulation(
        double fx, double fy,     /// input: the focal lengths (same for both cameras)
        double cx, double cy,     /// input: the principal point (same for both cameras)
        double s,                 /// input: the skew factor (same for both cameras)
        const std::vector<Vector2D> &points_0,  /// input: 2D image points in the 1st image.
        const std::vector<Vector2D> &points_1,  /// input: 2D image points in the 2nd image.
        std::vector<Vector3D> &points_3d,       /// output: reconstructed 3D points
        Matrix33 &R,   /// output: 3 by 3 matrix, which is the recovered rotation of the 2nd camera
        Vector3D &t    /// output: 3D vector, which is the recovered translation of the 2nd camera
    ) const {

    std::cout << "\nTODO: I am going to implement the triangulation() function in the following file:" << std::endl
              << "\t    - triangulation_method.cpp\n\n";

    //--------------------------------------------------------------------------------------------------------------

    /// @validity_check of the input

    // Check if the focal lengths are positive
    if (fx <=0 || fy <= 0) {
        std::cout << ">> The focal lengths should be positive!" << std::endl;
        return false;
    }
    std::cout << ">> The focal lengths are positive. Proceeding..." << std::endl;

    // Check if the principal point is within the image
    if (cx < 0 || cx >= fx || cy < 0 || cy >= fy) {
        std::cout << ">> The principal point should be within the image!" << std::endl;
        return false;
    }
    std::cout << ">> The principal point is within the image. Proceeding..." << std::endl;

    /** Check if the format of the 2D image points is correct
     * This check is performed to see if the number of coordinates per point is correct
     * It also accepts input for 2D image points with 3 coordinates if the third coordinate is 1
     */
    for (int i = 0; i < points_0.size(); i++) {
        if (points_0[i].size() > 2 || points_1[i].size() > 2) {
            std::cout << ">> The format of the 2D image points should be correct!" << std::endl;
            return false;
        }
    }
    std::cout << ">> The format of the 2D image points is correct. Proceeding..." << std::endl;

    // Check if the 2D image points of both images are equal in size and have at least 8 points
    std::cout << ">> length of points_0 = " << points_0.size() << std::endl;

    if (points_0.size() != points_1.size() && points_0.size() < 8) {
        std::cout << ">> The 2D image points of both images should be equal in size and have at least 8 points!"
                  << std::endl;
        return false;
    }
    else {
        std::cout << ">> The two 2D image points are equal in size and have at least 8 points. Proceeding..." << std::endl;
    }

    // Initialize the number of points variable
    int num_of_points = points_0.size();

    /**
     * Prepare W matrix
     *
     * W is an N x 9 matrix derived from N >= 8 correspondences
     * (In practice, it is often better to use more than eight correspondences and create a
        larger W matrix because it reduces the affects of noisy measurements)
     *
     */

    Matrix W(num_of_points, 9, 0.0); // Create a matrix of size num_of_points x 9 with all elements as 0.0

    for (int i = 0; i < num_of_points; i++) {
        double u = points_0[i][0];
        double up = points_1[i][0];
        double v = points_0[i][1];
        double vp = points_1[i][1];

        W.set_row(i, {u * up, v * up, up, u * vp, v * vp, vp, u, v, 1});
    }

//    std::cout << "W: " << W << std::endl;

    /**
     * Minimize using SVD
     * The SVD decomposition of W is used to find the estimated fundamental matrix F^
     * The last column of V is the solution to the minimization problem
     */
    Vector f = minimize_using_svd(W);

    Matrix33 F_hat(f[0], f[1], f[2],
                   f[3], f[4], f[5],
                   f[6], f[7], f[8]);

    Matrix U_f(3, 3, 0.0);
    Matrix D_f(3, 3, 0.0);
    Matrix V_f(3, 3, 0.0);

    svd_decompose(F_hat, U_f, D_f, V_f);

    /// Enforce the @rank-2 constraint on F^ to obtain the final fundamental matrix F
    D_f.set_column(D_f.cols() - 1, Vector3D(0.0, 0.0, 0.0));

    Matrix F = U_f * D_f * V_f.transpose();

    /// Compute the essential matrix E;
    Matrix33 K(fx, s, cx,
               0, fy, cy,
               0, 0, 1);

    Matrix33 E = K.transpose() * F * K;

    /// Recover rotation R and t.
    Matrix U_e(3, 3, 0.0);
    Matrix D_e(3, 3, 0.0);
    Matrix V_e(3, 3, 0.0);

    svd_decompose(E, U_e, D_e, V_e);

    Matrix33 W3(0.0, -1.0, 0.0,
                1.0, 0.0, 0.0,
                0.0, 0.0, 1.0);

    /**
     * The four possible solutions for the relative pose are:
     */
    Matrix R1 = determinant(U_e * W3 * V_e.transpose()) * U_e * W3 * V_e.transpose();
    Matrix R2 = determinant(U_e * W3.transpose() * V_e.transpose()) * U_e * W3.transpose() * V_e.transpose();
    Vector3D t1 = U_e.get_column(U_e.cols() - 1);
    Vector3D t2 = -U_e.get_column(U_e.cols() - 1);

    // Print the results
    std::cout << ">> Possible solutions for the relative pose are:" << std::endl;
    std::cout << "\t" << ">> R1 = " << R1 << std::endl;
    std::cout << "\t" << ">> R2 = " << R2 << std::endl;

    std::cout << "\t" << ">> t1 = " << t1 << std::endl;
    std::cout << "\t" << ">> t2 = " << t2 << "\n" << std::endl;

    // TODO: Reconstruct 3D points. The main task is
    //      - triangulate a pair of image points (i.e., compute the 3D coordinates for each corresponding point pair)

    Matrix M_left(3, 4, 0.0);
    M_left[0][0] = 1.0;
    M_left[1][1] = 1.0;
    M_left[2][2] = 1.0;

    Matrix M1 = createM(K, R1, t1);
    Matrix M2 = createM(K, R1, t2);
    Matrix M3 = createM(K, R2, t1);
    Matrix M4 = createM(K, R2, t2);

    bool world1 = true;
    bool world2 = true;
    bool world3 = true;
    bool world4 = true;

    for (int i = 0; i < num_of_points; i++) {
        Matrix A_1 = createA(points_0[i], points_1[i], M1, M_left);
        Matrix A_2 = createA(points_0[i], points_1[i], M2, M_left);
        Matrix A_3 = createA(points_0[i], points_1[i], M3, M_left);
        Matrix A_4 = createA(points_0[i], points_1[i], M4, M_left);

        Vector3D P_1 = static_cast<Vector4D>(minimize_using_svd(A_1)).cartesian();
        Vector3D P_2 = static_cast<Vector4D>(minimize_using_svd(A_2)).cartesian();
        Vector3D P_3 = static_cast<Vector4D>(minimize_using_svd(A_3)).cartesian();
        Vector3D P_4 = static_cast<Vector4D>(minimize_using_svd(A_4)).cartesian();

        if (P_1[2] < 0) {
            world1 = false;
        }
        if (P_2[2] < 0) {
            world2 = false;
        }
        if (P_3[2] < 0) {
            world3 = false;
        }
        if (P_4[2] < 0) {
            world4 = false;
        }
    }

    std::cout << "world1 = " << world1 << std::endl;
    std::cout << "world2 = " << world2 << std::endl;
    std::cout << "world3 = " << world3 << std::endl;
    std::cout << "world4 = " << world4 << std::endl;

    for (int i = 0; i < num_of_points; i++) {

        if (world1) {
            R = R1;
            t = t1;
            Matrix A_final = createA(points_0[i], points_1[i], M1, M_left);
            Vector3D P_final = static_cast<Vector4D>(minimize_using_svd(A_final)).cartesian();
            points_3d.push_back(P_final);
        }
        if (world2) {
            R = R1;
            t = t2;
            Matrix A_final = createA(points_0[i], points_1[i], M2, M_left);
            Vector3D P_final = static_cast<Vector4D>(minimize_using_svd(A_final)).cartesian();
            points_3d.push_back(P_final);
        }
        if (world3) {
            R = R2;
            t = t1;
            Matrix A_final = createA(points_0[i], points_1[i], M3, M_left);
//            std::cout << "matrix A final " << A_final << std::endl;
            Vector3D P_final = static_cast<Vector4D>(minimize_using_svd(A_final)).cartesian();
            points_3d.push_back(P_final);
        }
        if (world4) {
            R = R2;
            t = t2;
            Matrix A_final = createA(points_0[i], points_1[i], M4, M_left);
            Vector3D P_final = static_cast<Vector4D>(minimize_using_svd(A_final)).cartesian();
            points_3d.push_back(P_final);
        }

    }


    // TODO: Don't forget to
    //          - write your recovered 3D points into 'points_3d' (so the viewer can visualize the 3D points for you);
    //          - write the recovered relative pose into R and t (the view will be updated as seen from the 2nd camera,
    //            which can help you check if R and t are correct).
    //       You must return either 'true' or 'false' to indicate whether the triangulation was successful (so the
    //       viewer will be notified to visualize the 3D points and update the view).
    //       There are a few cases you should return 'false' instead, for example:
    //          - function not implemented yet;
    //          - input not valid (e.g., not enough points, point numbers don't match);
    //          - encountered failure in any step.
    return points_3d.size() > 0;
}
