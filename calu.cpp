#pragma comment(linker, "/STACK:102400000,102400000")
#include "mex.h"
#include <iostream>
#include <algorithm>
#include <vector>
#include <cmath>
#include <cstring>
#include <cstdlib>
using namespace std;

const int maxs = 80;
const int maxn = 30000;
struct node {
	int VarNumber, CaseLength;
	int VarSample[maxn][maxs] = {};
	int VarRange[maxs][maxs] = {}, VarRangeLength[maxs] = {};
};
node LGObj;
int dag[maxs][maxs] = {}, vis[maxn] = {}, flag[maxs] = {}, cnt[maxs] = {}, mm, n;
int sample[maxs] = {}, flagg[maxs] = {};
double theta[maxs][65536][5] = {}, m[maxs][65536][5] = {}, ss[maxs][65536][5] = {}, datt[maxn][maxs] = {};
vector<int>father[maxs], path;
vector<int>id;


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {         //主函数
	int p, q, d, t, nC = 0, nF = 0, row, maxparent = 0, curnode;
	double *nodescore, lastscore, localmax, localscore = 0.0, *Temp;
	n = (int)mxGetScalar(mxGetField(prhs[0], 0, "VarNumber"));          //节点数
	LGObj.VarNumber = n;
	LGObj.CaseLength = (int)(mxGetScalar(mxGetField(prhs[0], 0, "CaseLength")));     //样本数量
	Temp = mxGetPr(mxGetField(prhs[0], 0, "VarRange"));           //每个节点可以取的值 例如1结点可以取值1 2 3那么VarRange[1][0]=1  VarRange[1][1]=2  VarRange[1][2]=3
	row = mxGetN(mxGetField(prhs[0], 0, "VarRange"));            //可以取值的数组行数

	for (p = 0; p < n; p++)
		for (q = 0; q < row; q++)
			LGObj.VarRange[p][q] = (int)Temp[p + q * n];

	Temp = mxGetPr(mxGetField(prhs[0], 0, "VarRangeLength"));      //每个节点可以取值的个数例如结点1可以取1 2 3那么值VarRangeLength[1]=3
	for (p = 0; p < n; p++) LGObj.VarRangeLength[p] = (int)Temp[p];
	Temp = mxGetPr(mxGetField(prhs[0], 0, "VarSample"));        //样本值

	for (p = 0; p < LGObj.CaseLength; p++) {
		for (q = 0; q < n; q++) {
			LGObj.VarSample[p][q] = (int)Temp[p + q * LGObj.CaseLength];
		}
		id.push_back(p);
	}
	Temp = mxGetPr(prhs[1]);
	for (p = 0; p < n; p++)
		for (q = 0; q < n; q++) {
			dag[p][q] = (int)Temp[p + q * n];
		}
	for (int i = 0; i < n; i++)
		for (int j = 0; j < n; j++) {
			if (dag[j][i]) father[i].push_back(j);
		}
	Temp = mxGetPr(prhs[2]);
	int ld = 50;
	for (p = 0; p < ld; p++) {
		for (q = 0; q < n; q++) {
			datt[p][q] = (int)Temp[p + q * ld];
		}
	}
	int kk = 0;
	mm = LGObj.CaseLength;
	memset(vis, 0, sizeof(vis));
	memset(theta, 0, sizeof(theta));
	memset(sample, 0, sizeof(sample));
	memset(ss, 0, sizeof(ss));
	for (int i = 0; i < n; i++) cnt[i] = 1;
	for (int i = 0; i < LGObj.CaseLength; i++) {
		for (int j = 0; j < n; j++) {
			int tem1 = 0;
			if (father[j].size() > 0) sample[father[j].size() - 1] = 1;
			else sample[0] = 1;
			for (int k = father[j].size() - 2; k >= 0; k--) sample[k] = LGObj.VarRangeLength[father[j][k + 1]] * sample[k + 1];
			for (int k = 0; k < father[j].size(); k++) tem1 += sample[k] * (LGObj.VarSample[i][father[j][k]] - 1);
			if (tem1 > 0) tem1--;
			theta[j][tem1][LGObj.VarSample[i][j] - 1] ++;
		}
	}

	for (int i = 0; i < ld; i++) {
		for (int j = 0; j < n; j++) {
			int tem1 = 0;
			for (int k = 0; k < father[j].size(); k++) tem1 += sample[k] * (datt[i][father[j][k]] - 1);
			if (tem1 > 0) tem1--;
			ss[j][tem1][(int)datt[i][j] - 1] ++;
		}
	}
	for (int j = 0; j < n; j++) {
		int tem1 = 0;
		if (father[j].size() > 0) sample[father[j].size() - 1] = 1;
		for (int k = father[j].size() - 2; k >= 0; k--) sample[k] = LGObj.VarRangeLength[father[j][k + 1]] * sample[k + 1];
		if (father[j].size() > 0) cnt[j] = sample[0] * LGObj.VarRangeLength[father[j][0]];
	}
	for (int i = 0; i < n; i++)
		for (int j = 0; j < cnt[i]; j++)
			for (int k = 0; k < LGObj.VarRangeLength[i]; k++) {
				m[i][j][k] = theta[i][j][k];
				theta[i][j][k] /= LGObj.CaseLength;
			}
	double score = 0;
	for (int u = 0; u < n; u++)
		for (int j = 0; j < LGObj.VarRangeLength[u]; j++)
			for (int k = 0; k < cnt[u]; k++)if (theta[u][k][j] > 0) {
				printf("%lf %lf\n", ss[u][k][j] / 50.0, theta[u][k][j]);
				score += fabs(ss[u][k][j] / 50.0 - theta[u][k][j])*(-log(theta[u][k][j]));
			}
	plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
	Temp = mxGetPr(plhs[0]);
	Temp[0] = score;
}