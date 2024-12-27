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
            void calDens();               //�����ܶ�
            void calPress();                 //����ѹ��
            void calAccelerate();           //������ٶ�
            void reVelocityAndPosition();    //�����ٶȺ�λ��
            void checkBound();                //���߽�����
            void reBlock();                //���¿�ID

            ParticleSystem2d& mPs;       
            spline m;  //����Ȩ�ص�������������
        };
    }
}

#endif
