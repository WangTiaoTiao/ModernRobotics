#include "ur5.h"
#include "rclcpp/rclcpp.hpp"
#include "rclcpp_action/rclcpp_action.hpp"
#include "mr_interface_pkg/action/move_arm.hpp"
#include "sensor_msgs/msg/joint_state.hpp"

using namespace robot;

class ArmMoveServer : public rclcpp::Node
{
public:
    using action_t = mr_interface_pkg::action::MoveArm;
    using topic_t = sensor_msgs::msg::JointState;

    ArmMoveServer() : Node("actionS_move_twist_node")
    {   
        using namespace std::placeholders;
        
        robot_ = std::make_shared<UR5>(); 

        thetalist_ = Eigen::VectorXf::Zero(6);
        thetalist_ << -0.085, -0.9, 1.5, -2.12, -1.47, 0;

        js_msg_ = std::make_shared<topic_t>();
        js_msg_->name = {"shoulder_1_joint", "shoulder_2_joint",  "elbow_joint", 
                            "wrist_1_joint",    "wrist_2_joint", "wrist_3_joint"};
        js_msg_->position = {0, 0, 0, 0, 0, 0};
        
        js_pub_ = this->create_publisher<topic_t>("/joint_states", 10);

        server_ = rclcpp_action::create_server<action_t>(
                    this, "action_move_twist",
                    std::bind(&ArmMoveServer::handle_goal, this, _1, _2),
                    std::bind(&ArmMoveServer::handle_cancel, this, _1),
                    std::bind(&ArmMoveServer::handle_acceptable, this, _1));
        
        timer_ = this->create_wall_timer(10ms, std::bind(&ArmMoveServer::send_js, this));
    
        RCLCPP_INFO(this->get_logger(), "Create actionS_move_twist_node !");
    }

    ~ArmMoveServer()
    {
        RCLCPP_INFO(this->get_logger(), "Destory actionS_move_twist_node !");
    }
private:

    rclcpp_action::GoalResponse handle_goal(const rclcpp_action::GoalUUID &uuid, std::shared_ptr<const action_t::Goal> goal)
    {   
        (void) uuid;

        RCLCPP_INFO(this->get_logger(), " Imcoming action request id =: %d", goal->id);

        Eigen::VectorXf goal_twist(6);
        for (int i = 0; i < 6; ++i) goal_twist(i) = goal->goal_twist[i];
        const Eigen::Matrix4f tf_goal = se3ToTF(VecTose3(goal_twist));

        thetalist_goal_ = Eigen::VectorXf::Zero(6);
        bool flag = robot_->NewtonIK(tf_goal, thetalist_goal_, 0.1, 0.1);

        if (flag) {
            for (auto t : thetalist_goal_) RCLCPP_INFO(this->get_logger(), " goal thetas : %f", t);
        } else {
            RCLCPP_ERROR(this->get_logger(), " Cannot Inverse");
            return rclcpp_action::GoalResponse::REJECT;
        }

        return rclcpp_action::GoalResponse::ACCEPT_AND_EXECUTE; 
    }

    rclcpp_action::CancelResponse handle_cancel(const std::shared_ptr<rclcpp_action::ServerGoalHandle<action_t>> goal_handle)
    {
        (void) goal_handle;
        RCLCPP_INFO(this->get_logger(), " Cancel Command Accept !");
        return rclcpp_action::CancelResponse::ACCEPT;
    }

    void handle_acceptable(const std::shared_ptr<rclcpp_action::ServerGoalHandle<action_t>> goal_handle)
    {   
        std::thread(std::bind(&ArmMoveServer::execute, this, goal_handle)).detach();
    }

    void execute(const std::shared_ptr<rclcpp_action::ServerGoalHandle<action_t>> goal_handle)
    {
        const float total_time = goal_handle->get_goal()->total_time;

        const Eigen::VectorXf thetalist_start = thetalist_;

        auto feedback_msg = std::make_shared<action_t::Feedback>();
        auto result_msg = std::make_shared<action_t::Result>();

        const int hz = 30;
        const float step = 1.0f / (float)hz;
        rclcpp::Rate rate(hz);
        for (float i = 0.0f; (i <= total_time) && rclcpp::ok(); i = i + step)
        {
            if (goal_handle->is_canceling())
            {
                goal_handle->canceled(result_msg);
                return;
            }

            thetalist_ = TrapezoidalTrajectory(thetalist_start, thetalist_goal_, 0.3, total_time, i, "vt");
            
            feedback_msg->cur_time = i; 
            feedback_msg->diff = (thetalist_ - thetalist_goal_).norm();

            goal_handle->publish_feedback(feedback_msg);

            rate.sleep();
        }

        if (rclcpp::ok())
        {
            for (int i = 0; i < 6; ++i) result_msg->final_thetalist.push_back(thetalist_[i]); 
            goal_handle->succeed(result_msg);
        }
    }
    
    void send_js() const {
        for (int j = 0; j < 6; ++j) js_msg_->position[j] = deg2rad(rad2deg(thetalist_[j]));
        js_msg_->header.stamp = this->now();
        js_pub_->publish(*js_msg_);
    }

    std::shared_ptr<UR5> robot_;
    rclcpp_action::Server<action_t>::SharedPtr server_;
    rclcpp::Publisher<topic_t>::SharedPtr js_pub_;
    Eigen::VectorXf thetalist_;
    Eigen::VectorXf thetalist_goal_;
    rclcpp::TimerBase::SharedPtr timer_;
    std::shared_ptr<topic_t> js_msg_;

};

int main(int argc, char ** argv)
{
    rclcpp::init(argc, argv);
    rclcpp::spin(std::make_shared<ArmMoveServer>());
    rclcpp::shutdown();
    return 0;
}