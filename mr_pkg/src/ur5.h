#pragma once

#include "mr_h.h"

namespace robot  
{

class UR5
{
public:

    UR5() {
        initialize();
        cout << "[UR5] Arm is ready !" << endl;
    }

    Matrix4f ForwardK(const Eigen::VectorXf &thetalist) const {
        return FKinSpace(M_, Slist_, thetalist);
    }

    MatrixXf JacobianS(const Eigen::VectorXf &thetalist) const {
        return JacobianSpace(Slist_, thetalist);
    }

    MatrixXf JacobianB(const Eigen::VectorXf &thetalist) const {
        return JacobianBody(Blist_, thetalist);
    }

    bool NewtonIK(const Matrix4f &t_sd, VectorXf &thetalist_0, const float &e_omg, const float &e_v) {
        bool flag = IKinSpace(Slist_, M_, t_sd, thetalist_0, e_omg, e_v);
        return flag;
    }

    float ConditionNumber(const Eigen::VectorXf &thetalist) const {
        Eigen::MatrixXf js = JacobianS(thetalist);
        Eigen::MatrixXf A = js * js.transpose();
        
        // u1 closer to 1 is better
        // Eigen::EigenSolver<Eigen::MatrixXf> solver(A);
        // Eigen::VectorXf eigen_values = solver.eigenvalues().real();
        // float u1 = eigen_values.maxCoeff() / eigen_values.minCoeff();

        // u3 larger is better
        float u3 = sqrt(A.determinant());
        return u3;
    }
    
    // VectorXf NumericalIK(const Matrix4f &t_sd) {

    //     float a2 = 0.425, a3 = 0.39225, d6=0.0823, d4 = 0.10915;
    //     Matrix3f rot = t_sd.block<3, 3>(0, 0);
    //     Vector4f p6{0, 0, -d6, 1};
    //     Vector3f p = (t_sd * p6).head<3>();
    //     float d = (p.squaredNorm()  - a2 * a2 - a3 * a3) / (2 * a2 * a3);
        
    //     VectorXf thetalist(6);
    //     thetalist(0) = atan2(p(1), p(0));
    //     thetalist(2) = atan2(sqrt(1 - d * d), d);
    //     thetalist(1) = atan2(p(2), sqrt(p.head<2>().norm())) - atan2(a3 * sin(thetalist(2)), a2 + a3 * cos(thetalist(2)));
    //     thetalist(3) = atan2(rot(2, 0), rot(0, 0));
    //     thetalist(4) = acos(() / d6)
    //     thetalist(5) = atan2(rot(1, 2), rot(1, 0));

    //     return thetalist;
    // }

    ~UR5() {}

private:

    void initialize() 
    {

        float l1 = 0.425, l2 = 0.39225, w1 = 0.13585, w2 = 0.093, w3 = 0.0823, h1 = 0.089159, h2 = 0.09465;

        Slist_.reserve(6);
        Eigen::VectorXf s(6);
        s << 0, 0, 1, 0, 0, 0;
        Slist_.emplace_back(s);
        s << 0, 1, 0, -h1, 0, 0;
        Slist_.emplace_back(s);
        s << 0, 1, 0, -h1, 0, l1;
        Slist_.emplace_back(s);
        s << 0, 1, 0, -h1, 0, l1 + l2;
        Slist_.emplace_back(s);
        s << 0, 0, -1, -w1, l1 + l2, 0;
        Slist_.emplace_back(s);
        s << 0, 1, 0, h2 - h1, 0, l1 + l2;
        Slist_.emplace_back(s);

        Mlist_.reserve(7);
        Mlist_.emplace_back(Eigen::Matrix4f{{1, 0, 0,  0}, {0, 1, 0,  0}, { 0,  0, 1,  h1}, {0, 0, 0, 1}});
        Mlist_.emplace_back(Eigen::Matrix4f{{0, 0, 1,  0}, {0, 1, 0, w1}, {-1,  0, 0,   0}, {0, 0, 0, 1}});
        Mlist_.emplace_back(Eigen::Matrix4f{{1, 0, 0,  0}, {0, 1, 0,  0}, { 0,  0, 1,  l1}, {0, 0, 0, 1}});
        Mlist_.emplace_back(Eigen::Matrix4f{{0, 0, 1,  0}, {0, 1, 0,  0}, {-1,  0, 0,  l2}, {0, 0, 0, 1}});
        Mlist_.emplace_back(Eigen::Matrix4f{{1, 0, 0,  0}, {0, 1, 0, w2}, { 0,  0, 1,   0}, {0, 0, 0, 1}});
        Mlist_.emplace_back(Eigen::Matrix4f{{1, 0, 0,  0}, {0, 1, 0,  0}, { 0,  0, 1,  h2}, {0, 0, 0, 1}});
        Mlist_.emplace_back(Eigen::Matrix4f{{1, 0, 0,  0}, {0, 0, 1, w3}, { 0, -1, 0,   0}, {0, 0, 0, 1}});

        M_ = Eigen::Matrix4f{{-1, 0, 0, l1 + l2}, {0, 0, 1, w1 + w2}, {0, 1, 0, h1 - h2}, {0, 0, 0, 1}};
        
        auto ad_tbs =  getAdT(TransInv(M_));
        for (auto S : Slist_) {
            Blist_.emplace_back(ad_tbs * S);
        }

        Glist_.reserve(6);
        Glist_.emplace_back(getG_cylinder(3.7,    0.06, 0.15));
        Glist_.emplace_back(getG_cylinder(8.393,  0.06, 0.56));
        Glist_.emplace_back(getG_cylinder(2.275,  0.06, 0.50));
        Glist_.emplace_back(getG_cylinder(1.219,  0.6,  0.12));
        Glist_.emplace_back(getG_cylinder(1.219,  0.6,  0.12));
        Glist_.emplace_back(getG_cylinder(0.1879, 0.6,  0.12));
    }

    Eigen::MatrixXf getG_cylinder(const float &m, const float &r, const float &l) const {
        Eigen::MatrixXf g = Eigen::MatrixXf::Zero(6, 6);
        float ixy = m * (3 * r * r + l * l) / 12;
        float iz = m * r * r * 0.5;
        g.diagonal() << ixy, ixy, iz, m, m, m;
        return g;
    }

    std::vector<Eigen::VectorXf> Slist_;
    std::vector<Eigen::VectorXf> Blist_;
    std::vector<Eigen::Matrix4f> Mlist_;
    std::vector<Eigen::MatrixXf> Glist_;
    Eigen::Matrix4f M_;

};

    


}