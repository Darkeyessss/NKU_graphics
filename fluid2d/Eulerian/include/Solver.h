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

            // ��װ����
            void solve();
            void advect();
            void project();

            void advVel();//�ٶȶ���
            void applyForce();// ʩ����
            void advTemp();//�¶ȶ���
            void advDens();//�ܶȶ���

            //����ϵͳ����
            void buildPremaA();
            void buildb(unsigned int numCells);
            void buildPrecon();

        protected:
            MACGrid2d& mGrid;

            MACGrid2d target;//Ŀ�����񣬴洢������һ����״̬

            compressed_matrix<double> A;//ϡ����󣬱�ʾ���Է������ϵ������
            vector<double> b;//���Է�������Ҷ���
            vector<double> precon;//Ԥ�������������ٵ��������

            //�����
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

            // ǰ�������������ʱ����
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

            // ��������������м�����
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

            // Ӧ��Ԥ������
            applyprecon2d(PremaA, precon, currentResidual, searchDirect);

            double residualNorm = inner_prod(searchDirect, currentResidual);

            for (int iteration = 0; iteration < maximumIterations; ++iteration) {
                vector<double> matrixVector = prod(PremaA, searchDirect);
                double step = residualNorm / inner_prod(searchDirect, matrixVector);

                // ����ѹ����
                for (unsigned int i = 0; i < preSolu.size(); ++i) {
                    preSolu(i) += step * searchDirect(i);
                }

                // ���²в�
                for (unsigned int i = 0; i < currentResidual.size(); ++i) {
                    currentResidual(i) -= step * matrixVector(i);
                }

                double currentNorm = norm_2(currentResidual);
                if (currentNorm < tolerance) {
                    return true; // �����ɹ�
                }

                // Ӧ��Ԥ���������µĲв�
                applyprecon2d(PremaA, precon, currentResidual, currentIntermediate);

                double newResidualNorm = inner_prod(currentIntermediate, currentResidual);
                double scaling = newResidualNorm / residualNorm;

                // ������������
                for (unsigned int i = 0; i < searchDirect.size(); ++i) {
                    searchDirect(i) = currentIntermediate(i) + scaling * searchDirect(i);
                }

                residualNorm = newResidualNorm;
            }

            return false; // δ��������������������
        }

    }
}

#endif //!__EULER_SOLVER_H__