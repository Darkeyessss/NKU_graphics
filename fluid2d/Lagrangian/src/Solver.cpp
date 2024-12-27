#include "Lagrangian/include/Solver.h"
#include "Global.h"
#include <iostream>
#include <algorithm>
using namespace std;

namespace FluidSimulation
{
    namespace Lagrangian2d
    {

        Solver::Solver(ParticleSystem2d& ps) : mPs(ps), m(ps.supportRadius){}

        void Solver::solve()
        {
            // TODO
            // Solves the fluid simulation by performing some steps, which may include:
            // 1. compute density 
            // 2. compute press
            // 3. compute accleration
            // 4. update velocity and position
            // 5. check boundary
            // 6. update block id
            // ...
            calDens();// 1. compute density 
            calPress();// 2. compute press
            calAccelerate();// 3. compute accleration
            reVelocityAndPosition();// 4. update velocity and position
            checkBound();// 5. check boundary
            reBlock();// 6. update block id
        }

        // �����ܶ�
        void Solver::calDens()
        {
            for (auto& particle : mPs.particles)
            {
                particle.density = 0.0f;

                for (const auto& offset : mPs.blockIdOffs)
                {
                    int neighborBlockId = particle.blockId + offset;
                    if (neighborBlockId < 0 || neighborBlockId >= static_cast<int>(mPs.blockExtens.size()))
                        continue;

                    const auto& block = mPs.blockExtens[neighborBlockId];
                    for (int j = block.x; j < block.y; ++j)
                    {
                        if (&particle == &mPs.particles[j])
                            continue;

                        glm::vec2 displacement = particle.position - mPs.particles[j].position;
                        float dist = glm::length(displacement);
                        if (dist > Lagrangian2dPara::supportRadius)
                            continue;

                        particle.density += m.Value(dist);
                    }
                }

                // ���������ܶ�
                particle.density *= (mPs.particleVolume * Lagrangian2dPara::density);
                particle.density = max<float>(particle.density, Lagrangian2dPara::density);
            }
        }

        //����ѹ��
        void Solver::calPress()
        {
            for (int i = 0; i < mPs.particles.size(); i++)
            {
                mPs.particles[i].pressure = Lagrangian2dPara::stiffness * (powf(mPs.particles[i].density / Lagrangian2dPara::density, Lagrangian2dPara::exponent) - 1.0f);
                mPs.particles[i].pressDivDens2 = mPs.particles[i].pressure / (mPs.particles[i].density * mPs.particles[i].density);
            }
        }

        // ������ٶ�
        void Solver::calAccelerate()
        {
            float viscosityCoefficient = 2.0f * (1.0f + 2.0f) * Lagrangian2dPara::viscosity;

            for (auto& particle : mPs.particles)
            {
                particle.accleration = glm::vec2(-Lagrangian2dPara::gravityX, -Lagrangian2dPara::gravityY);

                glm::vec2 cumulativeViscosityForce(0.0f);
                glm::vec2 cumulativePressureForce(0.0f);

                for (const auto& offset : mPs.blockIdOffs)
                {
                    int adjacentBlockId = particle.blockId + offset;

                    if (adjacentBlockId < 0 || adjacentBlockId >= static_cast<int>(mPs.blockExtens.size()))
                        continue;

                    const auto& neighborBlock = mPs.blockExtens[adjacentBlockId];

                    for (int neighborIdx = neighborBlock.x; neighborIdx < neighborBlock.y; ++neighborIdx)
                    {
                        if (&particle == &mPs.particles[neighborIdx])
                            continue;

                        glm::vec2 relativePosition = particle.position - mPs.particles[neighborIdx].position;
                        float distance = glm::length(relativePosition);

                        if (distance <= Lagrangian2dPara::supportRadius)
                        {
                            glm::vec2 velocityDifference = particle.velocity - mPs.particles[neighborIdx].velocity;
                            float velocityProjection = glm::dot(velocityDifference, relativePosition);

                            float denominator = (distance * distance) + (0.01f * Lagrangian2dPara::supportRadius * Lagrangian2dPara::supportRadius);

                            glm::vec2 weightGradient = m.Grad(relativePosition);

                            glm::vec2 viscosityTerm = (Lagrangian2dPara::density * mPs.particleVolume / mPs.particles[neighborIdx].density) * velocityProjection * weightGradient / denominator;
                            cumulativeViscosityForce += viscosityTerm;

                            glm::vec2 pressureTerm = mPs.particles[neighborIdx].density * (particle.pressDivDens2 + mPs.particles[neighborIdx].pressDivDens2) * weightGradient;
                            cumulativePressureForce += pressureTerm;
                        }
                    }
                }

                // ճ������ѹ���ӵ����ٶ�
                particle.accleration += cumulativeViscosityForce * viscosityCoefficient;
                particle.accleration -= cumulativePressureForce * mPs.particleVolume;
            }
        }


        //�����ٶȺ�λ��
        void Solver::reVelocityAndPosition()
        {
            for (int i = 0; i < mPs.particles.size(); i++)
            {
                // �����ٶ�
                mPs.particles[i].velocity += static_cast<float>(Lagrangian2dPara::dt) * mPs.particles[i].accleration;

                // �����ٶ�
                glm::vec2 clampedVelocity;
                for (int j = 0; j < 2; j++)
                {
                    clampedVelocity[j] = max(-Lagrangian2dPara::maxVelocity, min(mPs.particles[i].velocity[j], Lagrangian2dPara::maxVelocity));
                }
                mPs.particles[i].velocity = clampedVelocity;

                // ����λ��
                mPs.particles[i].position += static_cast<float>(Lagrangian2dPara::dt) * mPs.particles[i].velocity;
            }
        }

        //���߽�
        void Solver::checkBound()
        {
            for (int i = 0; i < mPs.particles.size(); i++)
            {
                bool invFlag = false;

                //�±߽�
                if (mPs.particles[i].position.x < mPs.lowerBound.x + Lagrangian2dPara::supportRadius)
                {
                    mPs.particles[i].velocity.x = abs(mPs.particles[i].velocity.x);
                    invFlag = true;
                }
                if (mPs.particles[i].position.y < mPs.lowerBound.y + Lagrangian2dPara::supportRadius)
                {
                    mPs.particles[i].velocity.y = abs(mPs.particles[i].velocity.y);
                    invFlag = true;
                }

                //�ϱ߽�
                if (mPs.particles[i].position.x > mPs.upperBound.x - Lagrangian2dPara::supportRadius)
                {
                    mPs.particles[i].velocity.x = -abs(mPs.particles[i].velocity.x);
                    invFlag = true;
                }
                if (mPs.particles[i].position.y > mPs.upperBound.y - Lagrangian2dPara::supportRadius)
                {
                    mPs.particles[i].velocity.y = -abs(mPs.particles[i].velocity.y);
                    invFlag = true;
                }

                // ����б߽練�����ٶ�˥��
                if (invFlag)
                {
                    mPs.particles[i].velocity *= Lagrangian2dPara::velocityAttenuation;
                }

                // ����λ�ú��ٶ��ڱ߽���
                glm::vec2 newPosition, newVelocity;
                for (int j = 0; j < 2; j++)
                {
                    newPosition[j] = max(mPs.lowerBound[j] + Lagrangian2dPara::supportRadius + Lagrangian2dPara::eps,
                        min(mPs.particles[i].position[j],
                            mPs.upperBound[j] - (Lagrangian2dPara::supportRadius + Lagrangian2dPara::eps)));
                    newVelocity[j] = max(-Lagrangian2dPara::maxVelocity,
                        min(mPs.particles[i].velocity[j], Lagrangian2dPara::maxVelocity));
                }
                mPs.particles[i].position = newPosition;
                mPs.particles[i].velocity = newVelocity;
            }
        }

        //���¿�ID
        void Solver::reBlock()
        {
            for (int i = 0; i < mPs.particles.size(); i++)
            {
                glm::vec2 deltaPos = mPs.particles[i].position - mPs.lowerBound;
                glm::vec2 blockPosition = glm::floor(deltaPos / mPs.blockSize);
                mPs.particles[i].blockId = static_cast<int>(blockPosition.y) * mPs.blockNum.x + static_cast<int>(blockPosition.x);
            }
        }
    }
}
