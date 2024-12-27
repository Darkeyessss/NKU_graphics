#include "Eulerian/include/Solver.h"
#include "Configure.h"
using namespace std;
using namespace glm;

namespace FluidSimulation
{
    namespace Eulerian2d
    {
        using namespace ublas;
        Solver::Solver(MACGrid2d& grid) : mGrid(grid)
        {
            mGrid.reset();      // 重置网格状态
            buildPremaA();      // 构建压力系数矩阵A
            buildPrecon();      // 构建预条件向量precon
        }

        void Solver::solve()
        {
            // TODO
            // Solves the fluid simulation by performing some steps, which may include:
            // 1. advection
            // 2. compute external forces
            // 3. projection
            // ...
            // 调用各个步骤的封装函数
            advect();           // 对流更新速度、温度和密度
            applyForce();       // 施加外部力
            project();          // 投影确保速度场的无散性
        }

        void Solver::advect()
        {
            target.reset();     // 重置目标网格
            advVel();           // 速度
            advTemp();          // 温度
            advDens();          // 密度
        }

        //确保速度场的无散性
        void Solver::project()
        {
            //解Ax = b获取压力p
            unsigned int numCells = Eulerian2dPara::theDim2d[0] * Eulerian2dPara::theDim2d[1];
            buildb(numCells);

            ublas::vector<double> p(numCells);

            psolve(A, precon, b, p, 500, 0.005);

            double scaleConstant = Eulerian2dPara::dt / Eulerian2dPara::airDensity;
            double pressureChange;

            FOR_EACH_LINE
            {
                if (mGrid.isValid(i, j, mGrid.X))
                {
                    if (mGrid.isSolidFace(i, j, mGrid.X))
                    {
                        target.mU(i, j) = 0.0;
                    }
                    else
                    {
                        int index1 = mGrid.getIndex(i, j);
                        int index2 = mGrid.getIndex(i - 1, j);
                        pressureChange = (p(index1) - p(index2)) / mGrid.cellSize;
                        double vel = mGrid.mU(i, j);
                        vel = vel - scaleConstant * pressureChange;
                        target.mU(i, j) = vel;
                    }
                }
                if (mGrid.isValid(i, j, mGrid.Y))
                {
                    // Hard-code boundary condition for now
                    if (mGrid.isSolidFace(i, j, mGrid.Y))
                    {
                        target.mV(i, j) = 0.0;
                    }
                    else
                    {
                        int index1 = mGrid.getIndex(i, j);
                        int index2 = mGrid.getIndex(i, j - 1);
                        pressureChange = (p(index1) - p(index2)) / mGrid.cellSize;
                        double vel = mGrid.mV(i, j);
                        vel = vel - scaleConstant * pressureChange;
                        target.mV(i, j) = vel;
                    }
                }
            }

            mGrid.mU = target.mU;
            mGrid.mV = target.mV;
            assert(mGrid.checkDivergence());
        }

        //求解压力
        void Solver::buildPremaA()
        {
            unsigned int numCells = mGrid.mSolid.data().size();
            A.resize(numCells, numCells, false);

            for (unsigned int row = 0; row < numCells; row++)
            {
                int ri, rj;
                mGrid.getCell(row, ri, rj); // Each row corresponds to a cell
                if (mGrid.isSolidCell(ri, rj))
                    continue;
                for (unsigned int col = 0; col < numCells; col++)
                {
                    int ci, cj;
                    mGrid.getCell(col, ci, cj); // Each col corresponds to a possible neighbor
                    if (mGrid.isSolidCell(ci, cj))
                        continue;
                    double coeff = mGrid.getPressureCoeffBetweenCells(ri, rj, ci, cj);
                    if (fabs(coeff) > 0.0001)
                    {
                        A(row, col) = coeff;
                    }
                }
            }
        }

#define VALA(r, c) (r!= -1 && c!= -1)? A(r, c) : 0
        // 构建预条件向量precon，加速共轭梯度求解器
        void Solver::buildPrecon()
        {
            precon.resize(A.size1());
            fill(precon.begin(), precon.end(), 0);

            double tau = 0.0; // Disable MIC(0) 0.97;
            for (unsigned int index = 0; index < A.size1(); index++)
            {
                int i, j;
                mGrid.getCell(index, i, j);
                if (mGrid.isSolidCell(i, j))
                    continue;

                int neighbori = mGrid.getIndex(i - 1, j);
                int neighborj = mGrid.getIndex(i, j - 1);
                double termi = neighbori != -1 ? A(index, neighbori) * precon(neighbori) : 0;
                double termj = neighborj != -1 ? A(index, neighborj) * precon(neighborj) : 0;

                double termii = 0;
                if (neighbori != -1)
                {
                    int neighborij = mGrid.getIndex(i - 1, j + 1);
                    double termii0 = (VALA(neighbori, neighborij));
                    double termii1 = precon(neighbori) * precon(neighbori);
                    termii = VALA(index, neighbori) * termii0 / termii1;
                }

                double termjj = 0;
                if (neighborj != -1)
                {
                    int neighborji = mGrid.getIndex(i + 1, j - 1);
                    double termjj0 = (VALA(neighborj, neighborji));
                    double termjj1 = precon(neighborj) * precon(neighborj);
                    termjj = VALA(index, neighborj) * termjj0 / termjj1;
                }

                double e = A(index, index) - termi * termi - termj * termj - tau * (termii + termjj);

                precon(index) = 1 / sqrt(e);
            }
        }

        //向量b，求解线性系统 Ax = b
        void Solver::buildb(unsigned int numCells)
        {
            b.resize(numCells);
            double constant = -(Eulerian2dPara::airDensity * mGrid.cellSize * mGrid.cellSize) / Eulerian2dPara::dt;
            for (unsigned int index = 0; index < numCells; index++)
            {
                int i, j;
                mGrid.getCell(index, i, j);
                if (!mGrid.isSolidCell(i, j))
                    b(index) = constant * mGrid.getDivergence(i, j);
                else
                    b(index) = 0;
            }
        }

        // 对流速度场，更新速度场的u和v分量
        void Solver::advVel()
        {
            FOR_EACH_LINE
            {
                // advect u
                if (mGrid.isValid(i, j, mGrid.X))
                {
                    vec2 pos = mGrid.getLeft(i, j);
                    vec2 newpos = mGrid.semiLagrangian(pos, Eulerian2dPara::dt);

                    vec2 newvel = mGrid.getVelocity(newpos);
                    target.mU(i, j) = newvel[mGrid.X];
                }
            // advect v
            if (mGrid.isValid(i, j, mGrid.Y))
            {
                vec2 pos = mGrid.getBottom(i, j);
                vec2 newpos = mGrid.semiLagrangian(pos, Eulerian2dPara::dt);
                vec2 newvel = mGrid.getVelocity(newpos);
                target.mV(i, j) = newvel[mGrid.Y];
            }
            }

            mGrid.mU = target.mU;
            mGrid.mV = target.mV;
        }

        void Solver::applyForce()
        {
            FOR_EACH_LINE
            {
                if (mGrid.isValid(i, j, mGrid.Y))
                {
                    vec2 pos = mGrid.getBottom(i, j);
                    double yforce = mGrid.getBoussinesqForce(pos);
                    double vel = mGrid.mV(i, j);
                    vel = vel + Eulerian2dPara::dt * yforce;
                    target.mV(i, j) = vel;
                }
            }
            mGrid.mV = target.mV;
        }
 
        // 对流温度场
        void Solver::advTemp()
        {
            FOR_EACH_CELL
            {
                vec2 pos = mGrid.getCenter(i, j);
                vec2 newpos = mGrid.semiLagrangian(pos, Eulerian2dPara::dt);
                double newt = mGrid.getTemperature(newpos);
                target.mT(i, j) = newt;
            }
            mGrid.mT = target.mT;
        }

        // 对流密度场
        void Solver::advDens()
        {
            FOR_EACH_CELL
            {
                vec2 pos = mGrid.getCenter(i, j);
                vec2 newpos = mGrid.semiLagrangian(pos, Eulerian2dPara::dt);
                double newd = mGrid.getDensity(newpos);
                target.mD(i, j) = newd;
            }
            mGrid.mD = target.mD;
        }

    }
}