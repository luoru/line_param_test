#include <iostream>
#include <Eigen/Core>
#include <Eigen/Dense>

using namespace std;

Eigen::Vector4d plucker2Cayley(const Eigen::Matrix<double, 6, 1> &plucker_line)
{
    Eigen::Vector3d v, n;
    v = plucker_line.tail(3);
    n = plucker_line.head(3);
    
    Eigen::Vector3d P0;
    double v_squred_norm;
    v_squred_norm = v.squaredNorm();
    P0 = v.cross(n) / v_squred_norm;

    Eigen::Vector3d vn = v / v_squred_norm;
    Eigen::Vector3d m;
    m = P0.cross(vn);
    double d = m.norm();

    Eigen::Matrix3d Q;
    Eigen::Vector3d vn_x_m;
    vn_x_m = vn.cross(m);
    Q.col(0) = vn;
    Q.col(1) = m / d;
    Q.col(2) = vn_x_m / vn_x_m.norm();

    Eigen::Matrix3d s_skew;
    Eigen::Matrix3d I = Eigen::Matrix3d::Identity();
    s_skew = (Q - I) * (Q + I).inverse();
    Eigen::Vector3d s;
    s(0) = s_skew(2, 1);
    s(1) = s_skew(0, 2);
    s(2) = s_skew(1, 0);

    Eigen::Vector4d cayley_param;
    cayley_param.tail(3) = s;
    cayley_param(0) = d;
    return cayley_param;
}

Eigen::Matrix3d skewVector(const Eigen::Vector3d &vec)
{
    Eigen::Matrix3d skew;
    skew << 0, -vec(2), vec(1),
            vec(2), 0, -vec(0),
            -vec(1), vec(0), 0;
    return skew;
}

Eigen::Matrix<double, 6, 1> cayley2Plucker(const Eigen::Vector4d &cayley_param)
{
    Eigen::Vector3d s = cayley_param.tail(3);
    double d = cayley_param(0);
    double s_squared_norm = s.squaredNorm();
    Eigen::Matrix3d Q = ((1 - s_squared_norm) * Eigen::Matrix3d::Identity() + 2 * skewVector(s) + 2 * s * s.transpose()) / (1 + s_squared_norm);
    Eigen::Vector3d vn = Q.col(0);
    Eigen::Vector3d m = d * Q.col(1);
    Eigen::Matrix<double, 6, 1> plucker_param;
    plucker_param.head(3) = m;
    plucker_param.tail(3) = vn;
    return plucker_param;
}

Eigen::Matrix<double, 6, 1> twoPoints2Plucker(const Eigen::Vector3d &pt1, const Eigen::Vector3d &pt2)
{
    Eigen::Vector3d v = pt2 - pt1;
    Eigen::Vector3d n = pt1.cross(pt2);
    double d = n.norm() / v.norm();
    Eigen::Matrix<double, 6, 1> plucker_param;
    plucker_param.tail(3) = v.normalized();
    plucker_param.head(3) = d * n.normalized();
    return plucker_param;
}


int main()
{
    Eigen::Vector3d pt1, pt2;
    pt1 << 1, 1, 0;
    pt2 << 1, 1, 5;
    Eigen::Matrix<double, 6, 1> plucker_param = twoPoints2Plucker(pt1, pt2);
    Eigen::Matrix<double, 4, 1> cayley_param = plucker2Cayley(plucker_param);
    Eigen::Matrix<double, 6, 1> plucker = cayley2Plucker(cayley_param);
    return 0;
}