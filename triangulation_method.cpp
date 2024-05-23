#include "triangulation.h"
#include "matrix_algo.h"
#include <easy3d/optimizer/optimizer_lm.h>


using namespace easy3d;

Vector minimize_using_svd(Matrix M) {
    int num_rows = M.rows();
    int num_cols = M.cols();
    Matrix U(num_rows, num_rows, 0.0);
    Matrix D(num_rows, num_cols, 0.0);
    Matrix V(num_cols, num_cols, 0.0);

    svd_decompose(M, U, D, V);
    return V.get_column(V.cols() - 1);
}

Matrix createM(Matrix K_i, Matrix R_i, Vector3D t_i) {
    Matrix M_i(3, 4, 0.0);
    M_i.set_column(0, R_i.get_column(0));
    M_i.set_column(1, R_i.get_column(1));
    M_i.set_column(2, R_i.get_column(2));
    M_i.set_column(3, t_i);
    Matrix result = K_i * M_i;
    return result;
}

Matrix createA(Vector2D p0, Vector2D p1, Matrix M_i, Matrix M_l) {
    Matrix A(4, 4, 0.0);
    A.set_row(0, p0[0] * M_l.transpose().get_column(2) - M_l.transpose().get_column(0));
    A.set_row(1, p0[1] * M_l.transpose().get_column(2) - M_l.transpose().get_column(1));
    A.set_row(2, p1[0] * M_i.transpose().get_column(2) - M_i.transpose().get_column(0));
    A.set_row(3, p1[1] * M_i.transpose().get_column(2) - M_i.transpose().get_column(1));
    return A;
}

/**
 * TODO: Finish this function for reconstructing 3D geometry from corresponding image points.
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
    /// NOTE: there might be multiple workflows for reconstructing 3D geometry from corresponding image points.
    ///       This assignment uses the commonly used one explained in our lecture.
    ///       It is advised to define a function for the sub-tasks. This way you have a clean and well-structured
    ///       implementation, which also makes testing and debugging easier. You can put your other functions above
    ///       triangulation(), or put them in one or multiple separate files.

    std::cout << "\nTODO: I am going to implement the triangulation() function in the following file:" << std::endl
              << "\t    - triangulation_method.cpp\n\n";

    std::cout << "[Liangliang]:\n"
                 "\tFeel free to use any provided data structures and functions. For your convenience, the\n"
                 "\tfollowing three files implement basic linear algebra data structures and operations:\n"
                 "\t    - Triangulation/matrix.h  Matrices of arbitrary dimensions and related functions.\n"
                 "\t    - Triangulation/vector.h  Vectors of arbitrary dimensions and related functions.\n"
                 "\t    - Triangulation/matrix_algo.h  Determinant, inverse, SVD, linear least-squares...\n"
                 "\tPlease refer to the above files for a complete list of useful functions and their usage.\n\n"
                 "\tIf you choose to implement the non-linear method for triangulation (optional task). Please\n"
                 "\trefer to 'Tutorial_NonlinearLeastSquares/main.cpp' for an example and some explanations.\n\n"
                 "\tIn your final submission, please\n"
                 "\t    - delete ALL unrelated test or debug code and avoid unnecessary output.\n"
                 "\t    - include all the source code (and please do NOT modify the structure of the directories).\n"
                 "\t    - do NOT include the 'build' directory (which contains the intermediate files in a build step).\n"
                 "\t    - make sure your code compiles and can reproduce your results without ANY modification.\n\n"
              << std::flush;


    //--------------------------------------------------------------------------------------------------------------
    // implementation starts ...

    // TODO: check if the input is valid (always good because you never known how others will call your function).
    std::cout << "length of points_0 = " << points_0.size() << std::endl;

    if (points_0.size() != points_1.size() && points_0.size() < 8) {
        std::cerr << "The 2D image points of both images should be equal in size and have at least 8 points!"
                  << std::endl;
    } else {
        std::cout << "The two 2D image points are equal in size and have at least 8 points. Proceeding..." << std::endl;
    }

    int num_of_points = points_0.size();

    // TODO: Estimate relative pose of two views. This can be subdivided into
    // prepare W matrix
    Matrix W(num_of_points, 9, 0.0);

    for (int i = 0; i < num_of_points; i++) {
        double u = points_0[i][0];
        double up = points_1[i][0];
        double v = points_0[i][1];
        double vp = points_1[i][1];

        W.set_row(i, {u * up, v * up, up, u * vp, v * vp, vp, u, v, 1});
    }

    //      - estimate the fundamental matrix F;
    Vector f = minimize_using_svd(W);

    Matrix33 F_hat(f[0], f[1], f[2],
                   f[3], f[4], f[5],
                   f[6], f[7], f[8]);

    Matrix U_f(3, 3, 0.0);
    Matrix D_f(3, 3, 0.0);
    Matrix V_f(3, 3, 0.0);

    svd_decompose(F_hat, U_f, D_f, V_f);

    D_f.set_column(D_f.cols() - 1, Vector3D(0.0, 0.0, 0.0));

    Matrix F = U_f * D_f * V_f.transpose();

    //      - compute the essential matrix E;
    Matrix33 K(fx, s, cx,
               0, fy, cy,
               0, 0, 1);

    Matrix33 E = K.transpose() * F * K;

    //      - recover rotation R and t.

    Matrix U_e(3, 3, 0.0);
    Matrix D_e(3, 3, 0.0);
    Matrix V_e(3, 3, 0.0);

    svd_decompose(E, U_e, D_e, V_e);

    Matrix33 W3(0.0, -1.0, 0.0,
                1.0, 0.0, 0.0,
                0.0, 0.0, 1.0);

    Matrix R1 = determinant(U_e * W3 * V_e.transpose()) * U_e * W3 * V_e.transpose();
    Matrix R2 = determinant(U_e * W3.transpose() * V_e.transpose()) * U_e * W3.transpose() * V_e.transpose();
    Vector3D t1 = U_e.get_column(U_e.cols() - 1);
    Vector3D t2 = -U_e.get_column(U_e.cols() - 1);

    std::cout << "R1 = " << R1 << std::endl;
    std::cout << "R2 = " << R2 << std::endl;
    std::cout << "t1 = " << t1 << std::endl;
    std::cout << "t2 = " << t2 << std::endl;

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


    for (int i = 0; i < num_of_points; i++) {
        // not too sure about changing points_0/points_1
        Matrix A_1 = createA(points_0[i], points_1[i], M1, M_left);
        Matrix A_2 = createA(points_0[i], points_1[i], M2, M_left);
        Matrix A_3 = createA(points_0[i], points_1[i], M3, M_left);
        Matrix A_4 = createA(points_0[i], points_1[i], M4, M_left);

        Vector3D P_1 = static_cast<Vector4D>(minimize_using_svd(A_1)).cartesian();
        Vector3D P_2 = static_cast<Vector4D>(minimize_using_svd(A_2)).cartesian();
        Vector3D P_3 = static_cast<Vector4D>(minimize_using_svd(A_3)).cartesian();
        Vector3D P_4 = static_cast<Vector4D>(minimize_using_svd(A_4)).cartesian();

//        Vector4D m1p1 = static_cast<Vector4D>(M1 * P_1).cartesian();
//        Vector4D m2p2 = M2 * P_2;
//        Vector4D m3p3 = M3 * P_3;
//        Vector4D m4p4 = M4 * P_4;

        if (i < 5) {
//            std::cout << "A1:\t\t" << A_1 << std::endl;
//            std::cout << "A2:\t\t" << A_2 << std::endl;
//            std::cout << "P1:\t\t" << P_1 << std::endl;
//            std::cout << "P2:\t\t" << P_2 << std::endl;
            std::cout << "p " << i << P_1 << std::endl;
            std::cout << "p " << i << P_2 << std::endl;
            std::cout << "p " << i << P_3 << std::endl;
            std::cout << "p " << i << P_4 << std::endl;
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
