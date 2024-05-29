#include "triangulation.h"
#include "matrix_algo.h"
#include <easy3d/optimizer/optimizer_lm.h>

using namespace easy3d;


Vector minimize_using_svd(Matrix M) {
    /** Minimise using SVD (from matrix_algo.h Easy3D) and return the last column of V
    * Description: This function takes a matrix M and performs SVD decomposition on it.
    * input: Matrix M
    * output: Vector - last column of V in carterian coordinates
    */
    int num_rows = M.rows();
    int num_cols = M.cols();

    try {
        Matrix U(num_rows, num_rows, 0.0);
        Matrix D(num_rows, num_cols, 0.0);
        Matrix V(num_cols, num_cols, 0.0);

        svd_decompose(M, U, D, V);

        // Print the matrix
    //    std::cout << "Matrix v: " << V.get_column(V.cols() - 1) << std::endl;
    //    std::cout << "-----------------" << std::endl;

        Vector V_last = V.get_column(V.cols() - 1);

        return V_last;

    } catch (...) {
        // If any error occurs, throw an exception
        throw std::runtime_error("Failed to minimise using SVD");
    }
}


Matrix createM(Matrix K_i, Matrix R_i, Vector3D t_i) {
    /** Create M
     * Description: This function creates the M matrix
     * Input: Matrix K_i - calibration matrix, Matrix R_i - rotation matrix, Vector3D t_i - translation vector
     * Output: Matrix M_i - projection matrix
     */
    Matrix M_i(3, 4, 0.0);

    try {
        M_i.set_column(0, R_i.get_column(0));
        M_i.set_column(1, R_i.get_column(1));
        M_i.set_column(2, R_i.get_column(2));

        M_i.set_column(3, t_i);

        // Apply the calibration matrix
        M_i = K_i * M_i;

        return M_i;

    } catch (...) {
        // If any error occurs, throw an exception
        throw std::runtime_error("Failed to create M matrix");
    }
}


Matrix createA(Vector2D p0, Vector2D p1, Matrix M_i, Matrix M_left) {
    /** Create A
    * Description: This function creates the A matrix
    * Input: Vector2D p0 - 2D image point in the first image, Vector2D p1 - 2D image point in the second image,
    *          Matrix M_i - projection matrix of the second camera, Matrix M_l - projection matrix of the first camera
    * Output: Matrix A - 4x4 matrix
    */
    // Inistialise A matrix 4x4 with 0.0
    Matrix A(4, 4, 0.0);

    try {
        A.set_row(0, p0[0] * M_left.transpose().get_column(2) - M_left.transpose().get_column(0));
        A.set_row(1, p0[1] * M_left.transpose().get_column(2) - M_left.transpose().get_column(1));
        A.set_row(2, p1[0] * M_i.transpose().get_column(2) - M_i.transpose().get_column(0));
        A.set_row(3, p1[1] * M_i.transpose().get_column(2) - M_i.transpose().get_column(1));

        return A;
    } catch (...) {
        throw std::runtime_error("Failed to create A matrix");
    }
}


Matrix createW(std::vector<Vector2D> p0_norms, std::vector<Vector2D> p1_norms) {
    int N = p0_norms.size();
    Matrix W(N, 9, 0.0); // Create a matrix of size num_of_points x 9 with all elements as 0.0

    for (int i = 0; i < N; i++) {
        double u = p0_norms[i][0];
        double up = p1_norms[i][0];
        double v = p0_norms[i][1];
        double vp = p1_norms[i][1];

        W.set_row(i, {u * up, v * up, up, u * vp, v * vp, vp, u, v, 1});
    }
    return W;
}


Vector2D find_mean_image_point(std::vector<Vector2D> points_vector) {
    int N = points_vector.size();
    Vector2D mean_point = {0.0, 0.0};

    for (int i = 0; i < N; i++) {
        mean_point += points_vector[i];
    }

    mean_point /= N;
    return mean_point;
}


double find_scaling_factor(std::vector<Vector2D> points) {
    int N = points.size();
    Vector2D mean_point = find_mean_image_point(points);

    double mean_distance = 0.0;

    for (int i = 0; i < N; i++) {
        mean_distance += (points[i] - mean_point).norm();
    }

    mean_distance /= N;

    /// Calculate the scaling factor
    double scale = sqrt(2) / mean_distance;
    return scale;
}


Matrix33 createT(std::vector<Vector2D> points) {
    // find mean point of image points
    Vector2D mean = find_mean_image_point(points);

    // find scaling factor
    double scale = find_scaling_factor(points);

    // construct T
    Matrix33 T(scale, 0, -scale * mean[0],
                 0, scale,  -scale * mean[1],
                 0, 0,     1);
    return T;
}


std::vector<Vector2D> normalize_points(std::vector<Vector2D> points, Matrix33 T) {
    /// Normalise the 2D image points
    int N = points.size();
    std::vector<Vector2D> points_normalised;

    for (int i = 0; i < N; i++) {

        Vector3D new_point = T * points[i].homogeneous();

//        std::cout << "new_point: " << new_point << std::endl;

        points_normalised.push_back(new_point.cartesian());
    }
    return points_normalised;
}


Matrix createF(Matrix W, Matrix T0, Matrix T1) {
    Vector f = minimize_using_svd(W);

    // recast f into 3x3 matrix
    Matrix33 F_hat(f[0], f[1], f[2],
                   f[3], f[4], f[5],
                   f[6], f[7], f[8]);

    Matrix U_f(3, 3, 0.0);
    Matrix D_f(3, 3, 0.0);
    Matrix V_f(3, 3, 0.0);

    svd_decompose(F_hat, U_f, D_f, V_f);

    /// Enforce the @rank2 constraint on F_hat to obtain the @fundamental_matrix F
    D_f.set_column(D_f.cols() - 1, Vector3D(0.0, 0.0, 0.0));

    /// Fundamental matrix F is responsible for the epipolar geometry between two views
    Matrix F = U_f * D_f * V_f.transpose();

    /// @De_normalise the fundamental matrix F
    Matrix F_result = T1.transpose() * F * T0;

    return F_result;
}


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

    //--------------------------------------------------------------------------------------------------------------

    /// @validity_check of the input

    // Check if the focal lengths are positive
    if (fx <=0 || fy <= 0) {
        std::cout << ">> The focal lengths should be positive!" << std::endl;
        return false;
    }
    std::cout << ">> The focal lengths are positive. \t\t\tProceeding..." << std::endl;

    // Check if the principal point is within the image
    if (cx < 0 || cx >= fx || cy < 0 || cy >= fy) {
        std::cout << ">> The principal point should be within the image!" << std::endl;
        return false;
    }
    std::cout << ">> The principal point is within the image. \t\t\tProceeding..." << std::endl;


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
    std::cout << ">> The format of the 2D image points is correct. \t\tProceeding..." << std::endl;

    // Check if the 2D image points of both images are equal in size and have at least 8 points
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

    /** @Normalise the values of the 2D image points
     * Translate them to the origin and scale them to have a mean distance of sqrt(2) from the origin
     */

    /// Create the normalisation matrices
    std::cout << "Creating T matrices" << std::endl;
    Matrix33 T_0 = createT(points_0);
    Matrix33 T_1 = createT(points_1);
    std::cout << "\t\t\t\t\t\tSucces! Proceeding..." << std::endl;

    /// Normalise the 2D image points
    std::cout << "Normalizing points" << std::endl;
    std::vector<Vector2D> points_0_normalised = normalize_points(points_0, T_0);
    std::vector<Vector2D> points_1_normalised = normalize_points(points_1, T_1);
    std::cout << "\t\t\t\t\t\tSucces! Proceeding..."<< std::endl;

    /// Prepare W matrix
    std::cout << "Creating W matrix" << std::endl;
    /**
     * W is an N x 9 matrix derived from N >= 8 correspondences
     * (In practice, it is often better to use more than eight correspondences and create a
        larger W matrix because it reduces the affects of noisy measurements)
     */
    Matrix W = createW(points_0_normalised, points_1_normalised);
    std::cout << "\t\tsize W: (" << W.rows() << "," << W.cols() << ")" << std::endl;
    std::cout << "\t\t\t\t\t\tSucces! Proceeding..."<< std::endl;

    /// Prepare F matrix
    std::cout << "Creating F matrix" << std::endl;
    /**
     * Minimise W using SVD to find the @estimated_fundamental_matrix_F_hat
     * The SVD decomposition of W is used to find the estimated fundamental matrix F^
     * The last column of V is the solution to the minimization problem
     *
     * @warning The calculated fundamental matrix and its approximation are for normalised points
     * and they require de-normalisation to obtain the final fundamental matrix
     */
    Matrix F = createF(W, T_0, T_1);
    std::cout << "\t\tsize F: (" << F.rows() << "," << F.cols() << ")" << std::endl;
    std::cout << "\t\t\t\t\t\tSucces! Proceeding..."<< std::endl;

    /// Prepare E matrix
    std::cout << "Creating E matrix" << std::endl;
    /** Compute @Essential_matrix E
     * The Essential matrix is a 3x3 matrix that encapsulates the geometric information between two calibrated views in stereo vision.
     * It has five degrees of freedom (3 for rotation, 2 for translation) and is defined up to scale.
     *
     * The rank-2 constraint comes from the fact that the Essential matrix is singular, i.e., its determinant is zero.
     * This is a result of the epipolar constraint in stereo vision, which states that all corresponding points lie on the same plane.
     */
    Matrix33 K(fx, s, cx,
               0, fy, cy,
               0, 0, 1);
    Matrix33 E = K.transpose() * F * K;
    std::cout << "\t\tsize E: (" << E.rows() << "," << E.cols() << ")" << std::endl;
    std::cout << "\t\t\t\t\t\tSucces! Proceeding..."<< std::endl;



    /// Recover rotation R and t.
    Matrix U_e(3, 3, 0.0);
    Matrix D_e(3, 3, 0.0);
    Matrix V_e(3, 3, 0.0);

    svd_decompose(E, U_e, D_e, V_e);

    Matrix33 W3(0.0, -1.0, 0.0,
                1.0, 0.0, 0.0,
                0.0, 0.0, 1.0);


    Matrix33 Z(0.0, 1.0, 0.0,
               -1.0, 0.0, 0.0,
               0.0, 0.0, 0.0);


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


    Matrix M_left(3, 4, 0.0);
    M_left[0][0] = 1.0;
    M_left[1][1] = 1.0;
    M_left[2][2] = 1.0;
    M_left = K * M_left;

    Matrix M1 = createM(K, R1, t1);
    Matrix M2 = createM(K, R1, t2);
    Matrix M3 = createM(K, R2, t1);
    Matrix M4 = createM(K, R2, t2);

    /// Initiate flags for 4 potential world coordinates
    std::vector<bool> world_flags(4, false);

    /// Initiate counters for the number of points in the world
    std::vector<int> world_counters(4,0);

    for (int i=0; i < num_of_points; i++){

        // Construct A matrix for each case
        Matrix A_1 = createA(points_0[i], points_1[i], M1, M_left);
        Matrix A_2 = createA(points_0[i], points_1[i], M2, M_left);
        Matrix A_3 = createA(points_0[i], points_1[i], M3, M_left);
        Matrix A_4 = createA(points_0[i], points_1[i], M4, M_left);

        // Construct the 3D points from camera a
        Vector3D P_1a = static_cast<Vector4D>(minimize_using_svd(A_1)).cartesian();
        Vector3D P_2a = static_cast<Vector4D>(minimize_using_svd(A_2)).cartesian();
        Vector3D P_3a = static_cast<Vector4D>(minimize_using_svd(A_3)).cartesian();
        Vector3D P_4a = static_cast<Vector4D>(minimize_using_svd(A_4)).cartesian();

        // construct the 3D points from camera b
        Vector3D P_1b = R1 * P_1a + t1;
        Vector3D P_2b = R1 * P_2a + t2;
        Vector3D P_3b = R2 * P_3a + t1;
        Vector3D P_4b = R2 * P_4a + t2;

        // Check if the z-coordinate of the 3D point is positive for both cameras
        if (P_1a[2] > 0 && P_1b[2] > 0){
            world_counters[0]++;
        }

        if (P_2a[2] > 0 && P_2b[2] > 0){
            world_counters[1]++;
        }

        if (P_3a[2] > 0 && P_3b[2] > 0){
            world_counters[2]++;
        }

        if (P_4a[2] > 0 && P_4b[2] > 0){
            world_counters[3]++;
        }
    }

    // Find the world version with the maximum count of positive z values
    int max_value = world_counters[0];
    int max_index = 0;

    for (int i = 1; i < world_counters.size(); i++) {
        if (world_counters[i] > max_value) {
            max_value = world_counters[i];
            max_index = i;
        }
    }

    // Set the world flags
    for (int i = 0; i < world_flags.size(); i++) {
        if (i == max_index) {
            world_flags[i] = true;
        }
    }

    for (int i = 0; i < world_counters.size(); i++) {
        std::cout << "world_counters[" << i << "] = " << world_counters[i] << std::endl;
    }

    for (int i = 0; i <world_flags.size(); i++){
        std::cout << "world_flags[" << i << "] = " << world_flags[i] << std::endl;
    }

    // Reconstruct the 3D points
    for (int i = 0; i < num_of_points; i++) {

        if (world_flags[0]) {
            R = R1;
            t = t1;
            Matrix A_final = createA(points_0[i], points_1[i], M1, M_left);
            Vector3D P_final = static_cast<Vector4D>(minimize_using_svd(A_final)).cartesian();
            points_3d.push_back(P_final);
        }
        if (world_flags[1]) {
            R = R1;
            t = t2;
            Matrix A_final = createA(points_0[i], points_1[i], M2, M_left);
            Vector3D P_final = static_cast<Vector4D>(minimize_using_svd(A_final)).cartesian();
            points_3d.push_back(P_final);
        }
        if (world_flags[2]) {
            R = R2;
            t = t1;
            Matrix A_final = createA(points_0[i], points_1[i], M3, M_left);
            Vector3D P_final = static_cast<Vector4D>(minimize_using_svd(A_final)).cartesian();
            points_3d.push_back(P_final);
        }
        if (world_flags[3]) {
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
    //          - encountered failure in any step.

    return points_3d.size() > 0;
}
