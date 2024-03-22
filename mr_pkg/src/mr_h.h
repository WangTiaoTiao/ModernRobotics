#pragma once

#include <iostream>
#include <vector>
#include <algorithm>
#include <Eigen/Dense>
#include <iomanip>

using namespace Eigen;
using namespace std;

#define pi 3.1415926f
#define ZERO_TOL 1e-5f
#define cout cout << setprecision(4)

Matrix3f getSkew(const Vector3f& w) {
    Matrix3f skew_m{{    0, -w[2],  w[1]}, 
                    { w[2],     0, -w[0]},
                    {-w[1],  w[0],    0}};
    return skew_m;
}

void regularM(Matrix4f &in) {
    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) {
            in(i, j) = (abs(in(i, j)) < ZERO_TOL) ? 0 : in(i, j);
        }
    }
}

void regularM(MatrixXf &in) {
    int row = in.rows();
    int col = in.cols();
    for (int i = 0; i < row; ++i) {
        for (int j = 0; j < col; ++j) {
            in(i, j) = (abs(in(i, j)) < ZERO_TOL) ? 0 : in(i, j);
        }
    }
}

void TransToRP(const Matrix4f &tf, Matrix3f &r, Vector3f &p) {
    r = tf.block<3, 3>(0, 0);
    p = tf.block<3, 1>(0, 3);
}

bool isEqualM(const Matrix4f &M1, const Matrix4f &M2) {
    Matrix4f diff = M1 - M2;
    return (diff.norm() < 0.1) ? true : false;
}

bool isEqualM(const Matrix3f &r1, const Matrix3f &r2) {
    Matrix3f diff = r1 - r2;
    return (diff.norm() < 1e-2) ? true : false;
}

float cot(const float& theta) {
    return cos(theta)/sin(theta);
}

void getW(const Matrix3f& r, Vector3f& w, float &theta ) {
  if (r.isIdentity()) {
    cout << "Case: theta = zero, the solution is undefined" << endl;
    w = Vector3f::Identity();
  } else if (abs(r.trace()+1.0f) < 0.001f) {
    float coef = 1.0f/sqrt(2.0f * (1.0f + r(2,2)));
    cout << "Case: theta = pi, there are more than one possible solutions" << endl;
    // three fessiable equations and w or -w are both vaild solutions.
    w = Vector3f{r(0, 2) * coef, r(1, 2) * coef, (1 + r(2, 2)) * coef};
    theta = M_PI;
  } else {
    theta = acosf(0.5f * (r.trace() - 1.0f));
    Matrix3f w_mat = (r - r.transpose()) / (2.0f*sin(theta));
    w = Vector3f{w_mat(2, 1), w_mat(0, 2), w_mat(1, 0)};
  }
}

// Matrix3f getRotation(const Vector3f &w) {
//     Vector3f w_hat = w.normalized();
//     float theta = w.norm();
//     Matrix3f w_skew = getSkew(w_hat);
//     Matrix3f r = Matrix3f::Identity() + sin(theta) * w_skew +
//                 (1 - cos(theta)) * w_skew * w_skew;
//     return r;
// }

Matrix3f getRotation(const Vector3f &w_hat, const float& theta) {
    Matrix3f w_skew = getSkew(w_hat);
    Matrix3f r = Matrix3f::Identity() + sin(theta) * w_skew +
                (1 - cos(theta)) * w_skew * w_skew;
    return r;
}

void getScrew(const Matrix4f& tf, VectorXf& screw_aixs, float& theta) {
    Matrix3f r = tf.block<3, 3>(0, 0);
    Vector3f p = tf.block<3, 1>(0, 3);

    Vector3f w;
    Vector3f v;

    if (r == Matrix3f::Identity()) {
        w = Vector3f::Zero();
        v = p.normalized();
        theta = p.norm();
    } else {
        getW(r, w, theta);
        Matrix3f skew_w = getSkew(w);
        Matrix3f g_inv = (1/theta) * Matrix3f::Identity() - 0.5f * skew_w + 
                         (1/theta - 0.5f * cot(theta/2.0f)) * skew_w * skew_w;
        v = g_inv * p;
    } 
    screw_aixs << w, v;
    return;
};

void getTwist(const Matrix4f& tf, VectorXf& twist) {
    Matrix3f r = tf.block<3, 3>(0, 0);
    Vector3f p = tf.block<3, 1>(0, 3);
    float theta;
    Vector3f w;
    Vector3f v;

    if (r == Matrix3f::Identity()) {
        w = Vector3f::Zero();
        v = p.normalized();
        theta = p.norm();
    } else {
        getW(r, w, theta);
        Matrix3f skew_w = getSkew(w);
        Matrix3f g_inv = (1/theta) * Matrix3f::Identity() - 0.5f * skew_w + 
                         (1/theta - 0.5f * cot(theta/2.0f)) * skew_w * skew_w;
        v = g_inv * p;
    } 
    twist << w, v;
    twist = twist * theta;
    return;
};

Matrix4f getTF(const VectorXf& screw_aixs, const float& theta) {
    Vector3f w = screw_aixs.head<3>();
    Vector3f v = screw_aixs.segment<3>(3);
    Matrix4f tf = Matrix4f::Identity();
    if (w.isZero()) {
        v = v * theta;
    } else {
        tf.block<3,3>(0,0) = getRotation(w, theta);
        Matrix3f skew_w = getSkew(w);
        Matrix3f g = Matrix3f::Identity() * theta
                      + (1.0f - cos(theta)) * skew_w 
                      + (theta - sin(theta)) * skew_w * skew_w;
        v = g * v;
    }
    tf(0,3) = v(0);
    tf(1,3) = v(1);
    tf(2,3) = v(2);
    return tf;
}

MatrixXf getAdT(const Matrix4f& tf) {
    Matrix<float, 6, 6> adj = Matrix<float, 6, 6>::Zero();
    Matrix3f r = tf.block<3, 3>(0, 0);
    Vector3f p = tf.block<3, 1>(0, 3);
    adj.block<3,3>(0, 0) = r;
    adj.block<3,3>(3, 3) = r;
    adj.block<3,3>(3, 0) = getSkew(p) * r;
    return adj;
}

float deg2rad(const float& deg) {
    float rad = fmodf(deg * M_PI / 180.0f, 2 * M_PI);
    if (rad < -M_PI) {
        rad += 2 * M_PI; 
    } else if (rad > M_PI) {
        rad -= 2 * M_PI; // 将大于 pi 的值转换为对应的小于 pi 的值
    }
    return rad;
}

float rad2deg(const float& rad) {
    float deg = fmodf(rad * 180.0f / M_PI, 360.0f);
    if (deg < -180.0f) {
        deg += 360.0f; 
    } else if (deg > 180.0f) {
        deg -= 360.0f; // 将大于 pi 的值转换为对应的小于 pi 的值
    }
    return deg;
}

void coutV(const Vector3f& v) {
    cout << std::setprecision(3) << "<" << v[0] << "," << v[1] << "," << v[2] << ">" << endl;
}

Matrix4f FKinBody(const Matrix4f &M, const vector<VectorXf> &Blist, const vector<float> &thetalist) {
    Matrix4f ret = M;
    for (int i = 0; i < (int)thetalist.size(); ++i) {
        ret = ret * getTF(Blist[i], thetalist[i]);
    }
    regularM(ret);
    return ret;
}

Matrix4f FKinBody(const Matrix4f &M, const vector<VectorXf> &Blist, const VectorXf &thetalist) {
    Matrix4f ret = M;
    for (int i = 0; i < (int)thetalist.size(); ++i) {
        ret = ret * getTF(Blist[i], thetalist(i));
    }
    regularM(ret);
    return ret;
}

Matrix4f FKinSpace(const Matrix4f &M, const vector<VectorXf> &Slist, const vector<float> &thetalist) {
    Matrix4f ret = Matrix4f::Identity();
    for (int i = 0; i < (int)thetalist.size(); ++i) {
        ret = ret * getTF(Slist[i], thetalist[i]);
    }
    ret = ret * M;
    regularM(ret);
    return ret;
}

Matrix4f FKinSpace(const Matrix4f &M, const vector<VectorXf> &Slist, const VectorXf &thetalist) {
    Matrix4f ret = Matrix4f::Identity();
    for (int i = 0; i < thetalist.size(); ++i) {
        ret = ret * getTF(Slist[i], thetalist(i));
    }
    ret = ret * M;
    regularM(ret);
    return ret;
}

Vector3f getV(const Vector3f &q, const Vector3f &w, const float &h) {
    Vector3f v;
    v = (-w).cross(q) + h * w;
    return v;
}

void getScrewList(vector<VectorXf> &screwlist, const vector<Vector3f> &wlist, const vector<Vector3f> &qlist) {

    int n = wlist.size();
    VectorXf screw(6);
    Vector3f v;
    Vector3f w;
    for (int i = 0; i < n; ++i) {
        w = wlist[i];
        if (w.isZero()) {
            v = qlist[i];
        } else if (w.norm() != 1) {
            float h = w.norm();
            w = w.normalized();
            v = (-w).cross(qlist[i]) + w * h;
        } else {
            v = (-w).cross(qlist[i]);
        }
        screw << w, v;
        screwlist.push_back(screw);
    }
}

MatrixXf JacobianSpace(const vector<VectorXf> &ScrewList, const vector<float> &thetalist) {
    int n = ScrewList.size();
    MatrixXf jacob(6, n);

    vector<Matrix4f> tflist(n, Matrix4f::Identity());
    for (int i = 1; i < n; ++i) {
        tflist[i] = tflist[i-1] * getTF(ScrewList[i-1], thetalist[i-1]);
    }

    for (int i = 0; i < n; ++i) {
        auto adj = getAdT(tflist[i]);
        jacob.block<6, 1>(0, i) = adj * ScrewList[i];
    }
    return jacob;
}

MatrixXf JacobianSpace(const vector<VectorXf> &ScrewList, const VectorXf &thetalist) {
    int n = ScrewList.size();
    MatrixXf jacob(6, n);

    vector<Matrix4f> tflist(n, Matrix4f::Identity());
    for (int i = 1; i < n; ++i) {
        tflist[i] = tflist[i-1] * getTF(ScrewList[i-1], thetalist(i-1));
    }

    for (int i = 0; i < n; ++i) {
        auto adj = getAdT(tflist[i]);
        jacob.block<6, 1>(0, i) = adj * ScrewList[i];
    }
    return jacob;
}

MatrixXf JacobianBody(const vector<VectorXf> &ScrewList, const vector<float> &thetalist) {
    int n = ScrewList.size();
    MatrixXf jacob(6, n);

    vector<Matrix4f> tflist(n, Matrix4f::Identity());
    for (int i = n - 2; i >= 0; --i) {
        tflist[i] = tflist[i+1] * (getTF(ScrewList[i+1], thetalist[i+1])).inverse();
    }
    for (int i = 0; i < n; ++i) {
        auto adj = getAdT(tflist[i]);
        jacob.block<6, 1>(0, i) = adj * ScrewList[i];
    }
    regularM(jacob);
    return jacob;
}

MatrixXf JacobianBody(const vector<VectorXf> &ScrewList, const VectorXf &thetalist) {
    int n = ScrewList.size();
    MatrixXf jacob(6, n);

    vector<Matrix4f> tflist(n, Matrix4f::Identity());
    for (int i = n - 2; i >= 0; --i) {
        tflist[i] = tflist[i+1] * (getTF(ScrewList[i+1], thetalist(i+1))).inverse();
    }
    for (int i = 0; i < n; ++i) {
        auto adj = getAdT(tflist[i]);
        jacob.block<6, 1>(0, i) = adj * ScrewList[i];
    }
    regularM(jacob);
    return jacob;
}

Matrix4f TransInv(const Matrix4f &tf) {
    Matrix3f r = tf.block<3, 3>(0, 0);
    Vector3f p = tf.block<3, 1>(0, 3);

    Matrix4f inv_tf = Matrix4f::Identity();
    inv_tf.block<3, 3>(0, 0) = r.transpose();
    inv_tf.block<3, 1>(0, 3) = -r.transpose() * p;

    return inv_tf;
}

bool nearZero(const float &x) {
    return abs(x) < ZERO_TOL;
}

Matrix3f rot2so3(const Matrix3f &r) {
    Matrix3f so3 = Matrix3f::Zero();
    Vector3f omg = Vector3f::Zero();
    Vector3f x = Vector3f::Zero();
    if(isEqualM(r, Matrix3f::Identity())) {
        return so3;
    } else if (nearZero(r.trace() + 1.0f)) {
        if (!nearZero(1 + r(2, 2))) {
            x << r(0, 2), r(1, 2), 1 + r(2, 2);
            omg = (1 / sqrtf(2 * (1 + r(2, 2)))) * x;
        } else if (!nearZero(1 + r(1, 1))) {
            x << r(0, 1), 1 + r(1, 1), r(2, 1);
            omg = (1 / sqrtf(2 * (1 + r(1, 1)))) * x;
        } else {
            x << 1 + r(0, 0), r(1, 0), r(2, 0);
            omg = (1 / sqrtf(2 * (1 + r(0, 0)))) * x;
        }
        so3 = getRotation(omg, pi);
    } else {
        float theta = acosf((r.trace() - 1) / 2);
        so3 = theta * (1 / (2 * sin(theta))) * (r - r.transpose());
    }
    return so3;
}

Matrix4f tf2se3(const Matrix4f &tf) {
    Matrix3f r;
    Vector3f p;
    TransToRP(tf, r, p);

    Matrix4f se3 = Matrix4f::Zero();
    Matrix3f so3 = rot2so3(r);

    if (isEqualM(so3, Matrix3f::Zero())) {
        se3.block<3,1>(0,3) = p;
    } else {
        float theta = acosf((r.trace() - 1) / 2);
        se3.block<3,3>(0,0) = so3;
        se3.block<3,1>(0,3) = ((Matrix3f::Identity() - so3 / 2) + (1 / theta - cot(theta / 2) / 2) * so3 * so3 / theta) * p;
    }
    regularM(se3);
    return se3;
}

VectorXf se3ToVec(const Matrix4f &se3) {
    VectorXf twist(6);
    twist << se3(2, 1), se3(0, 2), se3(1, 0), se3(0, 3), se3(1, 3), se3(2, 3);
    return twist;
}

Matrix4f VecTose3(const VectorXf &twist) {
    Matrix4f se3 = Matrix4f::Identity();
    Vector3f omg = twist.head<3>();
    Vector3f v   = twist.segment<3>(3);

    se3.block<3,3>(0, 0) = getSkew(omg);
    se3.block<3,1>(0, 3) = v;
    return se3;
}

/*
% Takes a 3x3 skew-symmetric matrix (an element of so(3)).
% Returns the corresponding 3-vector (angular velocity).
*/
Vector3f so3ToVec(const Matrix3f &so3) {
    Vector3f omg;
    omg(0) = so3(2, 1);
    omg(1) = so3(0, 2);
    omg(2) = so3(1, 0);
    return omg;
}

/*
% Takes a 3x3 so(3) representation of exponential coordinates.
% Returns R in SO(3) that is achieved by rotating about omghat by theta 
% from an initial orientation R = I.
*/
Matrix3f so3ToRot(const Matrix3f &so3){
    Matrix3f rot = Matrix3f::Identity();
    Vector3f omg = so3ToVec(so3);
    float theta = omg.norm();
    if (!nearZero(theta)) {
        Matrix3f omg_mat = so3 / theta;
        rot = Matrix3f::Identity() + sin(theta) * omg_mat
            + (1 - cos(theta)) * omg_mat * omg_mat;
    }
    return rot;
}

/*
% Takes a se(3) representation of exponential coordinates.
% Returns a T matrix in SE(3) that is achieved by traveling along/about the 
% screw axis S for a distance theta from an initial configuration T = I.
*/
Matrix4f se3ToTF(const Matrix4f &se3) {
    Vector3f omg = so3ToVec(se3.block<3,3>(0,0));
    Matrix4f tf = Matrix4f::Identity();
    float theta = omg.norm();
    if (nearZero(theta)) {
        tf.block<3, 1>(0, 3) = se3.block<3, 1>(0, 3);
    } else {
        Matrix3f omg_mat = se3.block<3,3>(0,0) / theta;
        Vector3f v = se3.block<3,1>(0,3) / theta;
        tf.block<3,3>(0,0) = so3ToRot(se3.block<3,3>(0,0));
        Matrix3f g = Matrix3f::Identity() * theta
                   + (1.0f - cos(theta)) * omg_mat 
                   + (theta - sin(theta)) * omg_mat * omg_mat;
        tf.block<3, 1>(0, 3) = g * v;      
    }
    return tf;
}

bool IKinBody(const vector<VectorXf> &Blist, const Matrix4f &M, const Matrix4f &t_sd, VectorXf &thetalist_0, const float &e_omg, const float &e_v ) {
    int max_iter = 20, n = thetalist_0.size();
    VectorXf thetalist = thetalist_0;

    Matrix4f t_bs = TransInv(FKinBody(M, Blist, thetalist));
    Matrix4f t_bd = t_bs * t_sd;
    VectorXf twist_b = se3ToVec(tf2se3(t_bd));

    float cur_e_w = (twist_b.head<3>()).norm();
    float cur_e_v = (twist_b.segment<3>(3)).norm();

    for (int i = 0; i < max_iter; ++i) {
        if (cur_e_v <= e_v && cur_e_w <= e_omg) {
            thetalist_0 = thetalist;
            cout << "Inverse Found !" << endl;
            return true;
        }
        auto j_b = JacobianBody(Blist, thetalist);
        auto j_b_inv = j_b.completeOrthogonalDecomposition().pseudoInverse();
        thetalist = thetalist + j_b_inv * twist_b;
        for (int i = 0; i < n; ++i) thetalist[i] = deg2rad(rad2deg(thetalist[i]));

        t_bs = TransInv(FKinBody(M, Blist, thetalist));
        t_bd = t_bs * t_sd;
        twist_b = se3ToVec(tf2se3(t_bd));

        cur_e_w = (twist_b.head<3>()).norm();
        cur_e_v = (twist_b.segment<3>(3)).norm();
        cout << i << " : " << cur_e_w  << " - " << cur_e_v << endl;
    }
    cout << "[IKinBody] Inverse Not Found !" << endl;
    return false;
}

/*
% Takes Slist: The joint screw axes in the space frame when the manipulator
%              is at the home position, in the format of a matrix with the
%              screw axes as the columns,
%       M: The home configuration of the end-effector,
%       T: The desired end-effector configuration Tsd,
%       thetalist0: An initial guess of joint angles that are close to 
%                   satisfying Tsd,
%       eomg: A small positive tolerance on the end-effector orientation 
%             error. The returned joint angles must give an end-effector 
%             orientation error less than eomg,
%       ev: A small positive tolerance on the end-effector linear position 
%           error. The returned joint angles must give an end-effector 
%           position error less than ev.
% Returns thetalist: Joint angles that achieve T within the specified 
%                    tolerances,
%         success:
*/
bool IKinSpace(const vector<VectorXf> &Slist, const Matrix4f &M, const Matrix4f &t_sd, VectorXf &thetalist_0, const float &e_omg, const float &e_v ) {
    int max_iter = 20, n = thetalist_0.size();
    VectorXf thetalist = thetalist_0;

    Matrix4f t_sb = FKinSpace(M, Slist, thetalist);
    VectorXf twist_s = getAdT(t_sb) * se3ToVec(tf2se3(TransInv(t_sb) * t_sd));
    
    float cur_e_w = (twist_s.head<3>()).norm();
    float cur_e_v = (twist_s.segment<3>(3)).norm();

    for (int i = 0; i < max_iter; ++i) {
        if (cur_e_v <= e_v && cur_e_w <= e_omg) {
            thetalist_0 = thetalist;
            cout << "Inverse Found !" << endl;
            return true;
        }
        auto j_s = JacobianSpace(Slist, thetalist);
        auto j_s_inv = j_s.completeOrthogonalDecomposition().pseudoInverse();
        thetalist = thetalist + j_s_inv * twist_s;
        for (int i = 0; i < n; ++i) thetalist[i] = deg2rad(rad2deg(thetalist[i]));

        t_sb = FKinSpace(M, Slist, thetalist);
        twist_s = getAdT(t_sb) * se3ToVec(tf2se3(TransInv(t_sb) * t_sd));

        cur_e_w = (twist_s.head<3>()).norm();
        cur_e_v = (twist_s.segment<3>(3)).norm();
        cout << i << " : " << cur_e_w  << " - " << cur_e_v << endl;
    }
    cout << "[IKinBody] Inverse Not Found !" << endl;
    return false;
}

/*
% Takes V: 6-vector spatial velocity.
% Returns adV: The corresponding 6x6 matrix.
*/
MatrixXf getAdV(const VectorXf &V) {
    MatrixXf adv(6, 6);
    adv.setZero();

    Vector3f w = V.head<3>();
    Vector3f v = V.segment<3>(3);

    adv.block<3,3>(0,0) = getSkew(w);
    adv.block<3,3>(3,3) = getSkew(w);
    adv.block<3,3>(3,0) = getSkew(v); 

    regularM(adv);

    return adv;
}

/*
% Takes thetalist: n-vector of joint variables,
%       dthetalist: n-vector of joint rates,
%       ddthetalist: n-vector of joint accelerations,
%       g: Gravity vector g,
%       Ftip: Spatial force applied by the end-effector expressed in frame 
%             {n+1},
%       Mlist: List of link frames {i} relative to {i-1} at the home 
%              position,
%       Glist: Spatial inertia matrices Gi of the links,
%       Slist: Screw axes Si of the joints in a space frame, in the format
%              of a matrix with the screw axes as the columns.
% Returns taulist: The n-vector of required joint forces/torques.
*/
VectorXf InverseDynamic(const VectorXf &thetalist, const VectorXf &d_thetalist,   const VectorXf &dd_thetalist,  const Vector3f &g, 
                        const VectorXf &f_tip,     const vector<Matrix4f> &Mlist, const vector<MatrixXf> &Glist, const vector<VectorXf> &Slist)
{ 
    int n = thetalist.rows();
    Matrix4f Mi = Matrix4f::Identity();
    Matrix4f Ti = Matrix4f::Identity();
    vector<VectorXf> Ai(n);
    vector<VectorXf> Vi(n + 1, VectorXf(6));
    Vi[0] << 0, 0, 0, 0, 0, 0;
    vector<VectorXf> Vdi(n + 1, VectorXf(6));
    Vdi[0] << 0, 0, 0, -g;
    vector<MatrixXf> AdTi(n + 1, MatrixXf(6,6));
    AdTi[n] = getAdT(TransInv(Mlist[n]));
    VectorXf Fi = f_tip;
    VectorXf taulist(n);

    for (int i = 0; i < n; ++i) {

        // 8.50
        Mi = Mi * Mlist[i];
        Ai[i] = getAdT(TransInv(Mi)) * Slist[i];
        Ti = getTF(Ai[i], -thetalist[i]) * TransInv(Mlist[i]);

        // 8.51
        AdTi[i] = getAdT(Ti);
        Vi[i + 1] = AdTi[i] * Vi[i] + Ai[i] * d_thetalist[i];

        // 8.52
        Vdi[i + 1] = AdTi[i] * Vdi[i] + getAdV(Vi[i + 1]) * Ai[i] * d_thetalist[i] + Ai[i] * dd_thetalist[i];
    }

    for (int i = n - 1; i >= 0; --i) {

        // 8.53
        auto adj1 = AdTi[i + 1];
        auto adj2 = getAdV(Vi[i + 1]);
        Fi = adj1.transpose() * Fi + Glist[i] * Vdi[i + 1] - adj2.transpose() * (Glist[i] * Vi[i + 1]);
        
        // 8.54
        taulist[i] = Fi.transpose() * Ai[i];
    }
    return taulist;
}

/*
% Takes thetalist: A list of joint variables,
%       dthetalist: A list of joint rates,
%       Mlist: List of link frames i relative to i-1 at the home position,
%       Glist: Spatial inertia matrices Gi of the links,
%       Slist: Screw axes Si of the joints in a space frame, in the format
%              of a matrix with the screw axes as the columns,
% Returns c: The vector c(thetalist,dthetalist) of Coriolis and centripetal
%            terms for a given thetalist and dthetalist.
*/
VectorXf VelQuadraticForces(const VectorXf &thetalist, const VectorXf &d_thetalist,  const vector<Matrix4f> &Mlist, const vector<MatrixXf> &Glist, const vector<VectorXf> &Slist) {
    int n = thetalist.rows();
    VectorXf dd_thetalist(n);
    dd_thetalist.setZero();
    Vector3f g = Vector3f::Zero();
    VectorXf f_tip(6);
    f_tip.setZero();

    auto ret = InverseDynamic(thetalist, d_thetalist, dd_thetalist, g, f_tip, Mlist, Glist, Slist);
    return ret;
}

/*
% Takes thetalist: A list of joint variables,
%       g: 3-vector for gravitational acceleration,
%       Mlist: List of link frames i relative to i-1 at the home position,
%       Glist: Spatial inertia matrices Gi of the links,
%       Slist: Screw axes Si of the joints in a space frame, in the format
%              of a matrix with the screw axes as the columns.
% Returns grav: The joint forces/torques required to overcome gravity at 
%               thetalist
*/
VectorXf GravityForces(const VectorXf &thetalist, const Vector3f &g,  const vector<Matrix4f> &Mlist, const vector<MatrixXf> &Glist, const vector<VectorXf> &Slist) {
    int n = thetalist.rows();
    VectorXf d_thetalist(n);
    d_thetalist.setZero();
    VectorXf dd_thetalist(n);
    dd_thetalist.setZero();
    VectorXf f_tip(6);
    f_tip.setZero();

    auto ret = InverseDynamic(thetalist, d_thetalist, dd_thetalist, g, f_tip, Mlist, Glist, Slist);
    return ret;
}


/*
% Takes thetalist: A list of joint variables,
%       Ftip: Spatial force applied by the end-effector expressed in frame 
%             {n+1},
%       Mlist: List of link frames i relative to i-1 at the home position,
%       Glist: Spatial inertia matrices Gi of the links,
%       Slist: Screw axes Si of the joints in a space frame, in the format 
%              of a matrix with screw axes as the columns,
% Returns JTFtip: The joint forces and torques required only to create the 
%                 end-effector force Ftip.
*/
VectorXf EndEffectorForces(const VectorXf &thetalist, const VectorXf &f_tip,  const vector<Matrix4f> &Mlist, const vector<MatrixXf> &Glist, const vector<VectorXf> &Slist) {
    int n = thetalist.rows();
    VectorXf d_thetalist(n);
    d_thetalist.setZero();
    VectorXf dd_thetalist(n);
    dd_thetalist.setZero();
    Vector3f g = Vector3f::Zero();

    auto ret = InverseDynamic(thetalist, d_thetalist, dd_thetalist, g, f_tip, Mlist, Glist, Slist);
    return ret;
}

/*
% Takes thetalist: A list of joint variables,
%       Mlist: List of link frames i relative to i-1 at the home position,
%       Glist: Spatial inertia matrices Gi of the links,
%       Slist: Screw axes Si of the joints in a space frame, in the format
%              of a matrix with the screw axes as the columns.
% Returns M: The numerical inertia matrix M(thetalist) of an n-joint serial
%            chain at the given configuration thetalist.
*/
MatrixXf MassMatrix(const VectorXf &thetalist, const vector<Matrix4f> &Mlist, const vector<MatrixXf> &Glist, const vector<VectorXf> &Slist) {
    int n = thetalist.rows();
    MatrixXf M(n,n);
    VectorXf dd_thetalist(n);
    VectorXf d_thetalist(n);
    d_thetalist.setZero();
    Vector3f g = Vector3f::Zero();
    VectorXf f_tip(6);
    f_tip.setZero();

    for(int i = 0; i < n; ++i) {
        dd_thetalist.setZero();
        dd_thetalist(i) = 1;
        M.block(0, i, n, 1) = InverseDynamic(thetalist, d_thetalist, dd_thetalist, g, f_tip, Mlist, Glist, Slist);
    }
    return M;
}

/*
% Takes thetalist: A list of joint variables,
%       dthetalist: A list of joint rates,
%       taulist: An n-vector of joint forces/torques,
%       g: Gravity vector g,
%       Ftip: Spatial force applied by the end-effector expressed in frame 
%             {n+1},
%       Mlist: List of link frames i relative to i-1 at the home position,
%       Glist: Spatial inertia matrices Gi of the links,
%       Slist: Screw axes Si of the joints in a space frame, in the format
%              of a matrix with the screw axes as the columns,
% Returns ddthetalist: The resulting joint accelerations.
*/
VectorXf ForwardDynamics(const VectorXf &thetalist, const VectorXf &d_thetalist,   const VectorXf &taulist,  const Vector3f &g, 
                         const VectorXf &f_tip,     const vector<Matrix4f> &Mlist, const vector<MatrixXf> &Glist, const vector<VectorXf> &Slist) {
    auto M = MassMatrix(thetalist, Mlist, Glist, Slist);
    auto cf = VelQuadraticForces(thetalist, d_thetalist, Mlist, Glist, Slist);
    auto gf = GravityForces(thetalist, g, Mlist, Glist, Slist);
    auto ef = EndEffectorForces(thetalist, f_tip, Mlist, Glist, Slist);

    VectorXf dd_thetalist = M.inverse() * (taulist - cf - gf - ef);
    return dd_thetalist;
}

/*
% Takes Tf: Total time of the motion in seconds from rest to rest,
%       t: The current time t satisfying 0 < t < Tf.
% Returns s: The path parameter s(t) corresponding to a third-order 
%            polynomial motion that begins and ends at zero velocity.
*/
float CubicTimeScaling(const float &total_time, const float &cur_time){
    float x = (cur_time / total_time);
    return 3.0f * x * x - 2.0f * x * x * x;
}

/*
% Takes Tf: Total time of the motion in seconds from rest to rest,
%       t: The current time t satisfying 0 < t < Tf.
% Returns s: The path parameter s(t) corresponding to a fifth-order
%            polynomial motion that begins and ends at zero velocity and 
%            zero acceleration.
*/
float QuinticTimeScaling(const float &total_time, const float &cur_time){
    float x = (cur_time / total_time);
    return 10.0f * pow(x, 3) - 15.0f * powf(x, 4) + 6.0f * pow(x, 5);
}

float TrapezoidalTimeScaling(const float &v, const float &a, const float &total_time, const float &cur_time) {
    if (cur_time >= 0 && cur_time <= (v / a)) {
        return 0.5f * a * cur_time * cur_time;
    } else if (cur_time > (v / a) && cur_time <= (total_time - v/a)) {
        return v * cur_time - (v * v) / (2 * a);
    } else if (cur_time > (total_time - v/a) && cur_time <= total_time) {
        return (2 * a * v * total_time - 2 * v * v - 
                a * a * (cur_time - total_time) * (cur_time - total_time)) / (2 * a);
    } else {
        return -1;
    }
}

float TrapezoidalVA(const float &v, const float &a, const float &cur_time) {
    if ((v * v / a) > 1) return -1;
    float total_time = (a + v * v) / (v * a);
    return TrapezoidalTimeScaling(v, a, total_time, cur_time);
}

float TrapezoidalVT(const float &v, const float &total_time, const float &cur_time) {
    float vt = v * total_time;
    if (vt <= 1 || vt > 2) return -1;
    float a = (v * v) / (vt - 1);
    return TrapezoidalTimeScaling(v, a, total_time, cur_time);
}

float TrapezoidalAT(const float &a, const float &total_time, const float &cur_time) {
    float at = a * total_time * total_time;
    if (at < 4) return -1;
    float v = 0.5 * (a * total_time - sqrt(a) * sqrt(at - 4));
    return TrapezoidalTimeScaling(v, a, total_time, cur_time);
}


/*
% Takes thetastart: The initial joint variables,
%       thetaend: The final joint variables,
%       Tf: Total time of the motion in seconds from rest to rest,
%       N: The number of points N > 1 (Start and stop) in the discrete 
%          representation of the trajectory,
%       method: The time-scaling method, where 3 indicates cubic 
%               (third-order polynomial) time scaling and 5 indicates 
%               quintic (fifth-order polynomial) time scaling.
% Returns traj
*/
void JointTrajectory(vector<VectorXf> &traj, const VectorXf &start, const VectorXf &end, const float &total_time, const int &n, const int method = 3) {
    float s, time_gap = total_time / (n - 1);
    VectorXf cur = start;
    traj.clear();

    for (int i = 0; i < n; ++i) {
        if (method == 3) {
            s = CubicTimeScaling(total_time, time_gap * i);
        } else {
            s = QuinticTimeScaling(total_time, time_gap * i);
        }
        cur = start + s * (end - start);
        traj.emplace_back(cur);
    }
    return;
}

VectorXf TrapezoidalTrajectory(const VectorXf &start, const VectorXf &end, const float &in1, const float &in2, const float &cur_time, const string &method = "vt") {
    float total_time = (method == "va") ? (in2 + in1 * in1) / (in1 * in2) : in2;
    float s;

    if (method == "vt") {
        s = TrapezoidalVT(in1, total_time, cur_time);
    } else if (method == "at") {
        s = TrapezoidalAT(in1, total_time, cur_time);
    } else {
        s = TrapezoidalVA(in1, in2, cur_time);
    }
    if (s != -1) return start + s * (end - start);
    cout << "[TrapezoidalTrajectory] Error" << endl; 
    return start;
}


/*
% Takes Xstart: The initial end-effector configuration,
%       Xend: The final end-effector configuration,
%       Tf: Total time of the motion in seconds from rest to rest,
%       N: The number of points N > 1 (Start and stop) in the discrete
%          representation of the trajectory,
%       method: The time-scaling method, where 3 indicates cubic
%               (third-order polynomial) time scaling and 5 indicates 
%               quintic (fifth-order polynomial) time scaling.
% Returns traj: 
*/
void ScrewTrajectory(vector<Matrix4f> &traj, const Matrix4f &start, const Matrix4f &end, const float &total_time, const int &n, const int method = 3) {
    Matrix4f cur = start;
    float s, time_gap = total_time / (n - 1);
    traj.clear();
    traj.reserve(n);

    for (int i = 0; i < n; ++i) {
        if (method == 3) {
            s = CubicTimeScaling(total_time, time_gap * i);
        } else {
            s = QuinticTimeScaling(total_time, time_gap * i);
        }
        cur = start * se3ToTF(s * (tf2se3(TransInv(start) * end)));
        traj.emplace_back(cur);
    }
    return;
}

/*
% Takes Xstart: The initial end-effector configuration,
%       Xend: The final end-effector configuration,
%       Tf: Total time of the motion in seconds from rest to rest,
%       N: The number of points N > 1 (Start and stop) in the discrete 
%          representation of the trajectory,
%       method: The time-scaling method, where 3 indicates cubic 
%               (third-order polynomial) time scaling and 5 indicates 
%               quintic (fifth-order polynomial) time scaling.
% Returns traj:
*/
void CartesianTrajectory(vector<Matrix4f> &traj, const Matrix4f &start, const Matrix4f &end, const float &total_time, const int &n, const int method = 3) {
    Matrix4f cur = start;
    float s, time_gap = total_time / (n - 1);
    Matrix3f rot_start, rot_end;
    Vector3f p_start, p_end;
    TransToRP(start, rot_start, p_start);
    TransToRP(end, rot_end, p_end);

    traj.clear();
    traj.reserve(n);

    for (int i = 0; i < n; ++i) {
        if (method == 3) {
            s = CubicTimeScaling(total_time, time_gap * i);
        } else {
            s = QuinticTimeScaling(total_time, time_gap * i);
        }
        cur.block<3, 3>(0, 0) = rot_start * so3ToRot(rot2so3(rot_start.transpose() * rot_end) * s);
        cur.block<3, 1>(0, 3) = p_start + s * (p_end - p_start);
        traj.emplace_back(cur);
    }
    return;
} 

