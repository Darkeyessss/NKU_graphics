#pragma once
#ifndef __LAGRANGIAN_2D_SOLVER_H__
#define __LAGRANGIAN_2D_SOLVER_H__

#include "ParticleSystem2d.h"
#include "Configure.h"
#include "L2d.h"

namespace FluidSimulation
{

    namespace Lagrangian2d
    {
        class Solver
        {
        public:
            Solver(ParticleSystem2d& ps);

            void solve();

        private:
            void calDens();               //计算密度
            void calPress();                 //计算压力
            void calAccelerate();           //计算加速度
            void reVelocityAndPosition();    //更新速度和位置
            void checkBound();                //检查边界条件
            void reBlock();                //更新块ID

            ParticleSystem2d& mPs;       
            spline m;  //计算权重的三次样条函数
        };
    }
}

#endif
