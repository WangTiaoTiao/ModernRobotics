#pragma once

#include <iostream>
#include <vector>
#include <algorithm>
#include <Eigen/Dense>
#include <iomanip>

using namespace Eigen;
using namespace std;

#define pi 3.1415926f
#define ZERO_TOL 1e-6f
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

void tf2rp(const Matrix4f &tf, Matrix3f &r, Vector3f &p) {
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
    cout << theta << endl;
    Matrix3f w_mat = (r - r.transpose()) / (2.0f*sin(theta));
    w = Vector3f{w_mat(2, 1), w_mat(0, 2), w_mat(1, 0)};
  }
}

Matrix3f getRotation(const Vector3f &w) {
    Vector3f w_hat = w.normalized();
    float theta = w.norm();
    Matrix3f w_skew = getSkew(w_hat);
    Matrix3f r = Matrix3f::Identity() + sin(theta) * w_skew +
                (1 - cos(theta)) * w_skew * w_skew;
    return r;
}

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
    cout << twist.transpose() << '\n' << theta << endl;
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

Matrix4f getTF(const VectorXf& twist) {
    VectorXf screw_axis(6);
    Vector3f w = twist.head<3>();
    Vector3f v = twist.segment<3>(2);
    float theta;
    if (w.isZero()) {
        theta = v.norm();
        v = v.normalized();
    } else {
        theta = w.norm();
        w = w.normalized();
        v = v / theta;
    }
    screw_axis << w, v;
    return getTF(screw_axis, theta);
    
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
    tf2rp(tf, r, p);

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

VectorXf se3toVec(const Matrix4f &se3) {
    VectorXf twist(6);
    twist << se3(2, 1), se3(0, 2), se3(1, 0), se3(0, 3), se3(1, 3), se3(2, 3);
    return twist;
}

bool IKinBody(const vector<VectorXf> &Blist, const Matrix4f &M, const Matrix4f &t_sd, VectorXf &thetalist_0, const float &e_omg, const float &e_v ) {
    int max_iter = 20, n = thetalist_0.size();
    VectorXf thetalist = thetalist_0;

    Matrix4f t_bs = TransInv(FKinBody(M, Blist, thetalist));
    Matrix4f t_bd = t_bs * t_sd;
    VectorXf twist_b = se3toVec(tf2se3(t_bd));

    float cur_e_w = (twist_b.head<3>()).norm();
    float cur_e_v = (twist_b.segment<3>(3)).norm();

    for (int i = 0; i < max_iter; ++i) {
        if (cur_e_v <= e_v && cur_e_w <= e_omg) {
            thetalist_0 = thetalist;
            return true;
        }
        auto j_b = JacobianBody(Blist, thetalist);
        auto j_b_inv = j_b.completeOrthogonalDecomposition().pseudoInverse();
        thetalist = thetalist + j_b_inv * twist_b;
        for (int i = 0; i < n; ++i) thetalist[i] = deg2rad(rad2deg(thetalist[i]));

        t_bs = TransInv(FKinBody(M, Blist, thetalist));
        t_bd = t_bs * t_sd;
        twist_b = se3toVec(tf2se3(t_bd));

        cur_e_w = (twist_b.head<3>()).norm();
        cur_e_v = (twist_b.segment<3>(3)).norm();
        cout << i << " : " << cur_e_w  << " - " << cur_e_v << endl;
    }
    return false;
}

bool IKinSpace(const vector<VectorXf> &Slist, const Matrix4f &M, const Matrix4f &t_sd, VectorXf &thetalist_0, const float &e_omg, const float &e_v ) {
    int max_iter = 20, n = thetalist_0.size();
    VectorXf thetalist = thetalist_0;

    Matrix4f t_sb = FKinSpace(M, Slist, thetalist);
    VectorXf twist_s = getAdT(t_sb) * se3toVec(tf2se3(TransInv(t_sb) * t_sd));
    
    float cur_e_w = (twist_s.head<3>()).norm();
    float cur_e_v = (twist_s.segment<3>(3)).norm();

    for (int i = 0; i < max_iter; ++i) {
        if (cur_e_v <= e_v && cur_e_w <= e_omg) {
            thetalist_0 = thetalist;
            return true;
        }
        auto j_s = JacobianSpace(Slist, thetalist);
        auto j_s_inv = j_s.completeOrthogonalDecomposition().pseudoInverse();
        thetalist = thetalist + j_s_inv * twist_s;
        for (int i = 0; i < n; ++i) thetalist[i] = deg2rad(rad2deg(thetalist[i]));

        t_sb = FKinSpace(M, Slist, thetalist);
        twist_s = getAdT(t_sb) * se3toVec(tf2se3(TransInv(t_sb) * t_sd));

        cur_e_w = (twist_s.head<3>()).norm();
        cur_e_v = (twist_s.segment<3>(3)).norm();
        cout << i << " : " << cur_e_w  << " - " << cur_e_v << endl;
    }
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
        auto adj1 = AdTi[i + 1].transpose();
        auto adj2 = getAdV(Vi[i + 1]).transpose();
        Fi = adj1 * Fi + Glist[i] * Vdi[i + 1] - adj2 * (Glist[i] * Vi[i + 1]);

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
% Takes thetalist: n-vector of initial joint variables,
%       dthetalist: n-vector of initial joint rates,
%       taumat: An N x n matrix of joint forces/torques, where each row is 
%               the joint effort at any time step,
%       g: Gravity vector g,
%       Ftipmat: An N x 6 matrix of spatial forces applied by the 
%                end-effector (If there are no tip forces, the user should 
%                input a zero and a zero matrix will be used),
%       Mlist: List of link frames {i} relative to {i-1} at the home
%              position,
%       Glist: Spatial inertia matrices Gi of the links,
%       Slist: Screw axes Si of the joints in a space frame, in the format
%              of a matrix with the screw axes as the columns,
%       dt: The timestep between consecutive joint forces/torques,
%       intRes: Integration resolution is the number of times integration
%               (Euler) takes places between each time step. Must be an 
%               integer value greater than or equal to 1.
% Returns thetamat: The N x n matrix of robot joint angles resulting from 
%                   the specified joint forces/torques,
%         dthetamat: The N x n matrix of robot joint velocities.
*/
void ForwardDynamicsTrajectory() {
    return;
}

