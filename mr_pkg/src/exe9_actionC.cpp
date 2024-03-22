#include "ur5.h"
#include "rclcpp/rclcpp.hpp"
#include "rclcpp_action/rclcpp_action.hpp"
#include "mr_interface_pkg/action/move_arm.hpp"
#include "sensor_msgs/msg/joint_state.hpp"

using namespace robot;

class ArmMoveClient : public rclcpp::Node
{
public:
    using action_t = mr_interface_pkg::action::MoveArm;

    ArmMoveClient() : Node("actionC_move_twist_node")
    {   
        client_ = rclcpp_action::create_client<action_t>(this, "action_move_twist");
        RCLCPP_INFO(this->get_logger(), "Create actionC_move_twist_node !");
    }

    void send_goal(const vector<float> &twist, const int &id, const float &total_time)
    {
        using namespace std::chrono_literals;
        using namespace std::placeholders;

        RCLCPP_INFO(this->get_logger(), "Connecting to action server ----");
        if (!client_->wait_for_action_server(3s))
        {
            RCLCPP_ERROR(this->get_logger(), "Failed to Connect to action server");
            return;
        }

        auto goal = action_t::Goal();
        auto options = rclcpp_action::Client<action_t>::SendGoalOptions();

        RCLCPP_INFO(this->get_logger(), "Action server Connected ! Sending mission id : %d", id);
        
        goal.total_time = total_time;
        goal.id = id;
        for (int i = 0; i < 6; ++i) goal.goal_twist.push_back(twist[i]); 

        options.feedback_callback = std::bind(&ArmMoveClient::feedback_cb, this, _1, _2);
        options.goal_response_callback = std::bind(&ArmMoveClient::goal_response_cb, this, _1);
        options.result_callback = std::bind(&ArmMoveClient::result_cb, this, _1);
        auto future = client_->async_send_goal(goal, options);
    }
    
    ~ArmMoveClient()
    {
        RCLCPP_INFO(this->get_logger(), "Destory actionC_move_twist_node !");
    }
    
private:

    // rclcpp_action::Client<ActionT>::GoalResponseCallback = std::function<void (rclcpp_action::ClientGoalHandle<ActionT>::SharedPtr)>
    void goal_response_cb(const rclcpp_action::ClientGoalHandle<action_t>::SharedPtr &goal_handle)
    {
        if (!goal_handle)
        {
            RCLCPP_ERROR(this->get_logger(), " Goal was rejected by server");
        } else {
            RCLCPP_INFO(this->get_logger(), " Goal accepted by server, waiting for result");
        }
    }

    // std::function<void (typename ClientGoalHandle<ActionT>::SharedPtr,const std::shared_ptr<const Feedback>)>;
    void feedback_cb(rclcpp_action::ClientGoalHandle<action_t>::SharedPtr goal_handle, 
                     const std::shared_ptr<const action_t::Feedback> feedback)
    {
        (void) goal_handle;
        float cur_time = feedback->cur_time;
        float diff = feedback->diff;

        RCLCPP_INFO(this->get_logger(), "Current time : %.2f , with remained distance : %.2f",cur_time, diff);
    }

    // ResultCallback = std::function<void (const WrappedResult & result)>;
    void result_cb(const rclcpp_action::ClientGoalHandle<action_t>::WrappedResult &result)
    {   
        switch (result.code)
        {
            case rclcpp_action::ResultCode::ABORTED:
                RCLCPP_ERROR(this->get_logger(), "Goal was aborted");
                return;
            case rclcpp_action::ResultCode::CANCELED:
                RCLCPP_ERROR(this->get_logger(), "Goal was canceled");
                return;
            case rclcpp_action::ResultCode::SUCCEEDED:
                RCLCPP_INFO(this->get_logger(), " Successed !");
                for (size_t i = 0; i < 6; ++i) {
                    RCLCPP_INFO(this->get_logger(), " current thetas: %ld  : %f", i, result.result->final_thetalist[i]);
                }
                break;
            default:
                RCLCPP_ERROR(this->get_logger(), "Unkown result code");
                return;
        }

    }

    rclcpp_action::Client<action_t>::SharedPtr client_;
};

int main(int argc, char const *argv[])
{

    rclcpp::init(argc, argv);
    auto node = std::make_shared<ArmMoveClient>();

    UR5 robot;
    VectorXf thetalist(6);
    thetalist << 1, 2, 3, 4, 5, 6;
    int id = 1;
    float total_time = 4;
    auto twist = se3ToVec(tf2se3(robot.ForwardK(thetalist)));
    vector<float> vec;
    for (int i = 0; i < 6; ++i) vec.push_back(twist[i]);
    node->send_goal(vec, id, total_time);
    rclcpp::spin(node);
    rclcpp::shutdown();
    return 0;
}