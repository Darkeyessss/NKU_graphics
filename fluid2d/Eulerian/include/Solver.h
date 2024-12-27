#pragma once
#ifndef __EULERIAN_2D_SOLVER_H__
#define __EULERIAN_2D_SOLVER_H__

#pragma warning(disable: 4244 4267 4996)
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <cmath>
#include <algorithm>

#include "Eulerian/include/MACGrid2d.h"
#include "Global.h"


namespace FluidSimulation {
    namespace Eulerian2d {
        using namespace boost::numeric;
        using namespace ublas;
        class Solver {
        public:
            Solver(MACGrid2d& grid);

            // 封装函数
            void solve();
            void advect();
            void project();

            void advVel();//速度对流
            void applyForce();// 施加力
            void advTemp();//温度对流
            void advDens();//密度对流

            //线性系统构建
            void buildPremaA();
            void buildb(unsigned int numCells);
            void buildPrecon();

        protected:
            MACGrid2d& mGrid;

            MACGrid2d target;//目标网格，存储计算下一步的状态

            compressed_matrix<double> A;//稀疏矩阵，表示线性方程组的系数矩阵
            vector<double> b;//线性方程组的右端项
            vector<double> precon;//预处理向量，加速迭代求解器

            //求解器
            bool psolve(const compressed_matrix<double>& A,
                const vector<double>& precon,
                const vector<double>& b,
                vector<double>& p,
                int max_iter, double tol);

            void applyprecon2d(const compressed_matrix<double>& A,
                const vector<double>& precon,
                const vector<double>& r,
                vector<double>& z);
        };

        // Implementation of the functions
        inline void Solver::applyprecon2d(const compressed_matrix<double>& PremaA,
            const vector<double>& precon,
            const vector<double>& residual,
            vector<double>& intermediate) {

            unsigned int totalCells = residual.size();
            vector<double> temporary(totalCells, 0.0);

            // 前向遍历，计算临时向量
            for (unsigned int cellIdx = 0; cellIdx < totalCells; ++cellIdx) {
                int x, y;
                mGrid.getCell(cellIdx, x, y);

                int leftNeighbor = mGrid.getIndex(x - 1, y);
                int bottomNeighbor = mGrid.getIndex(x, y - 1);
                double contributionLeft = (leftNeighbor != -1) ? PremaA(cellIdx, leftNeighbor) * precon(leftNeighbor) * temporary(leftNeighbor) : 0.0;
                double contributionBottom = (bottomNeighbor != -1) ? PremaA(cellIdx, bottomNeighbor) * precon(bottomNeighbor) * temporary(bottomNeighbor) : 0.0;

                double adjustment = residual(cellIdx) - contributionLeft - contributionBottom;
                temporary(cellIdx) = adjustment * precon(cellIdx);
            }

            // 后向遍历，更新中间向量
            for (int cellIdx = static_cast<int>(totalCells) - 1; cellIdx >= 0; --cellIdx) {
                int x, y;
                mGrid.getCell(cellIdx, x, y);

                int rightNeighbor = mGrid.getIndex(x + 1, y);
                int topNeighbor = mGrid.getIndex(x, y + 1);
                double contributionRight = (rightNeighbor != -1) ? PremaA(cellIdx, rightNeighbor) * precon(cellIdx) * intermediate(rightNeighbor) : 0.0;
                double contributionTop = (topNeighbor != -1) ? PremaA(cellIdx, topNeighbor) * precon(cellIdx) * intermediate(topNeighbor) : 0.0;

                double adjustment = temporary(cellIdx) - contributionRight - contributionTop;
                intermediate(cellIdx) = adjustment * precon(cellIdx);
            }
        }


        inline bool Solver::psolve(const compressed_matrix<double>& PremaA,
            const vector<double>& precon,
            const vector<double>& rightHand,
            vector<double>& preSolu,
            int maximumIterations, double tolerance) {

            std::fill(preSolu.begin(), preSolu.end(), 0.0);
            vector<double> currentResidual = rightHand;
            vector<double> currentIntermediate = rightHand;
            vector<double> searchDirect = rightHand;

            // 应用预条件器
            applyprecon2d(PremaA, precon, currentResidual, searchDirect);

            double residualNorm = inner_prod(searchDirect, currentResidual);

            for (int iteration = 0; iteration < maximumIterations; ++iteration) {
                vector<double> matrixVector = prod(PremaA, searchDirect);
                double step = residualNorm / inner_prod(searchDirect, matrixVector);

                // 更新压力解
                for (unsigned int i = 0; i < preSolu.size(); ++i) {
                    preSolu(i) += step * searchDirect(i);
                }

                // 更新残差
                for (unsigned int i = 0; i < currentResidual.size(); ++i) {
                    currentResidual(i) -= step * matrixVector(i);
                }

                double currentNorm = norm_2(currentResidual);
                if (currentNorm < tolerance) {
                    return true; // 收敛成功
                }

                // 应用预条件器到新的残差
                applyprecon2d(PremaA, precon, currentResidual, currentIntermediate);

                double newResidualNorm = inner_prod(currentIntermediate, currentResidual);
                double scaling = newResidualNorm / residualNorm;

                // 更新搜索方向
                for (unsigned int i = 0; i < searchDirect.size(); ++i) {
                    searchDirect(i) = currentIntermediate(i) + scaling * searchDirect(i);
                }

                residualNorm = newResidualNorm;
            }

            return false; // 未能在最大迭代次数内收敛
        }

    }
}

#endif //!__EULER_SOLVER_H__