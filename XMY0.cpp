#pragma comment(linker, "/STACK:102400000,102400000")
#include "mex.h"
#include <iostream>
#include <algorithm>
#include <vector>
#include <cmath>
#include <cstring>
#include <cstdlib>
using namespace std;




typedef pair<double, int>P;
typedef double dlf;
typedef long long int LL;
const dlf inf = -999999999.0;
const int maxn = 30000;
const int maxs = 40;
const int maxmm = 10;
const dlf alpha = 0.9;
const dlf eps = 1e-6;
const dlf eps2 = 0.001;
dlf probilaty1[maxs][maxmm] = {}, probilaty2[maxs][maxmm][maxs][maxmm] = {}, probilaty3[maxs][maxmm][maxs][maxmm] = {}, entropy[maxs] = {}, betentropy[maxs][maxs] = {};
dlf probilaty4[maxs][maxmm][maxs][maxmm][maxs][maxmm] = {}, probilaty5[maxs][maxmm][maxs][maxmm][maxs][maxmm] = {}, MMI[maxs] = {};
int vis[maxs][maxs] = {}, vis2[maxs][6] = {}, tem1, tem2, tem3, tem4, tem5, tem6[5] = {}, tem7[6] = {}, dag[maxs][maxs] = {}, flagk = 0;
double pnum1[maxs][maxmm] = {}, pnum2[maxs][maxmm][maxs][maxmm] = {}, pnum3[maxs][maxmm][maxs][maxmm][maxs][maxmm] = {};
dlf maxx;
int par[maxs] = {}, rankk[maxs] = {}, viss[maxs][maxs][maxs] = {}, cntx = 0, cnty = 0, n;
int DAG[maxs][maxs] = {}, path[10] = {}, path2[10] = {}, cnpa1 = 0, cnpa2 = 0, value[maxs] = {}, Order[maxs] = {}, into[maxs] = {}, dist[maxs][maxs] = {}, vis3[maxs] = {};
vector<int>fa[maxs];
P anss[maxs*maxs] = {}, ess[maxs*maxs] = {};

//并查集
int found(int x) { return x == par[x] ? x : par[x] = found(par[x]); } //查找是否处于同一个联通块
void unite(int x, int y) {                                            //合并联通块
	x = found(x), y = found(y);
	if (x == y) return;
	if (rankk[x] < rankk[y]) par[x] = y;
	else {
		par[y] = x;
		if (rankk[x] == rankk[y]) rankk[x]++;
	}
}

inline dlf fabss(dlf x) { return x < 0 ? -x : x; }



struct edge {                                   //边与互信息
	int u, v;
	dlf w;
	bool operator < (const edge& rhs)const {
		return w > rhs.w;
	}
}e[maxs*maxs];
vector<int>g[maxs];


struct node {
	int VarNumber, CaseLength;
	int VarSample[maxn][maxs] = {};
	int VarRange[maxs][maxs] = {}, VarRangeLength[maxs] = {};
};
node LGObj;


dlf CH(node LGObj, int X, int *PAX, int cnt) {                                      //CH评分函数
	int N = LGObj.CaseLength, viss[LGObj.CaseLength + 100] = {}, TotalNumber = 0;
	int nX = LGObj.VarRangeLength[X];
	int RangeX[nX * 2] = {};
	dlf Frequency[2 * nX], GFunValue = 0.0;
	memset(viss, 0, sizeof(viss));
	for (int p = 0; p < N; p++) {
		int d = 1;
		for (int q = 0; q < cnt; q++)if (LGObj.VarSample[p][PAX[q]] == -1) {
			d = 0;
			break;
		}
		if (d == 0) {
			viss[p] = 1;
			TotalNumber++;
		}
	}
	viss[N] == TotalNumber;
	for (int i = 0; i < LGObj.VarRangeLength[X]; i++) RangeX[i] = LGObj.VarRange[X][i];
	int ri = nX, d = 0;
	while (d < N) {
		memset(Frequency, 0.0, sizeof(Frequency));
		while (d < N && viss[d] == 1) d++;
		if (d >= N) break;
		int t1;
		for (t1 = 0; t1 < nX; t1++)if (RangeX[t1] == LGObj.VarSample[d][X]) break;
		Frequency[t1] = 1;
		viss[d] = 1;
		dlf ParentValue[2 * cnt] = {};
		for (int i = 0; i < cnt - 1; i++) ParentValue[i] = LGObj.VarSample[d][PAX[i + 1]];
		d++;
		if (d >= N) break;
		for (int k = d; k < N; k++)if (!viss[k]) {
			int flag = 0;
			for (int i = 0; i < cnt - 1; i++)if (ParentValue[i] != LGObj.VarSample[k][PAX[i + 1]]) { flag = 1; break; }
			if (flag) continue;
			t1 = 0;
			while (RangeX[t1] != LGObj.VarSample[k][X] && t1 < nX) t1++;
			Frequency[t1]++;
			viss[k] = 1;
		}
		dlf sum = 0.0;
		for (int i = 0; i < nX; i++) sum = sum + Frequency[i];
		for (int k = 0; k < ri; k++)if (Frequency[k] != 0)
			GFunValue = GFunValue + lgamma(Frequency[k] + 1);
		GFunValue = GFunValue + lgamma(ri) - lgamma(sum + ri);
	}
	return GFunValue;
}




int dfs2(int s, int u, int fa, int d, int dd) {                     //枚举环的所用情况并进行评分  s为起始位置  u为当前位置  fa为上一个节点  d为当前度  dd为最大度
	if (d == dd && u == s) {
		dlf temm = 0.0;
		int PAX[10] = {}, cntt = 0;
		if (!path[0]) PAX[cntt++] = path2[1];
		if (path[cnpa2 - 1]) PAX[cntt++] = path2[cnpa2 - 1];
		temm = CH(LGObj, path2[0], PAX, cntt);
		for (int i = 1; i < cnpa2 - 1; i++) {
			cntt = 0;
			if (!path[i]) PAX[cntt++] = path2[i + 1];
			if (path[i - 1]) PAX[cntt++] = path2[i - 1];
			temm += CH(LGObj, path2[i], PAX, cntt);

		}
		cntt = 0;
		if (!path[cnpa2 - 1]) PAX[cntt++] = s;
		if (path[cnpa2 - 2]) PAX[cntt++] = path2[cnpa2 - 2];
		temm += CH(LGObj, u, PAX, cntt);
		if (temm > maxx) {
			maxx = temm; cnty = 1;
			for (int i = 0; i < cnpa2; i++) { tem6[i] = path[i]; tem7[i] = path2[i]; }
		}
		return true;
	}
	if (d >= dd)  return false;
	int flag = 0;
	path2[cnpa2++] = u; cnpa1++;              //path2记录遍历的节点编号  path1记录边的方向  0为u->v  1为v->u
	if (d == d - 1) {
		path[cnpa1 - 1] = 0;
		flag = max(flag, dfs2(s, s, u, d + 1, dd));    //若原来已经有方向则不改变原来规定的方向
		path[cnpa1 - 1] = 1;
		flag = max(flag, dfs2(s, s, u, d + 1, dd));
		cnpa1--; cnpa2--;
		return flag;
	}
	for (int v = 0; v < n; v++)if (v != u && v != fa) {
		path[cnpa1 - 1] = 0;
		if (d == dd - 2 && (!DAG[v][s] && !DAG[s][v])) continue;
		if (DAG[v][u]) flag = max(flag, dfs2(s, v, u, d + 1, dd));
		path[cnpa1 - 1] = 1;
		if (DAG[u][v]) flag = max(flag, dfs2(s, v, u, d + 1, dd));
	}
	cnpa1--; cnpa2--;
	return flag;
}



bool dfs3(int u) {
	vis2[u][0] = -1;
	for (int v = 0; v < n; v++)if (vis2[v][0] != 1 && DAG[u][v]) {     //判断是否存在环 vis2[u][d]为在度数为d的时候u节点的状态 -1代表处于栈中 1代表已经搜索完并且出栈0代表还未搜索
		if (vis2[v][0] == -1) return true;
		else if (dfs3(v)) return true;
	}
	vis2[u][0] = 1;
	return false;
}


bool calp(int *PAX, int x, int y, int cnt2, int d) {                           //判断隐藏父节点是服从 P(B|A,S)=P(B|C,A,S)
	if (d == cnt2) {
		int sum1 = 0, sum2 = 0, sum3 = 0, sum4 = 0;
		for (int i = 0; i < LGObj.CaseLength; i++) {
			int flagg = 0;
			for (int j = 0; j < cnt2; j++)if (value[j] != LGObj.VarSample[i][PAX[j]]) {
				flagg = 1; break;
			}
			if (!flagg) sum1++;
			flagg = 0;
			for (int j = 1; j < cnt2; j++)if (value[j] != LGObj.VarSample[i][PAX[j]]) {
				flagg = 1; break;
			}
			if (!flagg) sum2++;
			flagg = 0;
			for (int j = 2; j < cnt2; j++)if (value[j] != LGObj.VarSample[i][PAX[j]]) {
				flagg = 1; break;
			}
			if (!flagg) sum3++;
			flagg = 0;
			for (int j = 0; j < cnt2; j++)if (j != 1 && value[j] != LGObj.VarSample[i][PAX[j]]) {
				flagg = 1; break;
			}
			if (!flagg) sum4++;
		}
		if ((LL)sum1 * sum3 != (LL)sum2 * sum4) return false;
		return true;
	}
	for (int i = 0; i < LGObj.VarRangeLength[PAX[d]]; i++) {
		value[d] = LGObj.VarRange[PAX[d]][i];                   //value记录当前节点所取的值
		if (!calp(PAX, x, y, cnt2, d + 1)) return false;
	}
	return true;
}



void get_into_degree() {              //计算所有节点的入度数
	for (int j = 0; j < n; ++j) {
		into[j] = 0;
		for (int i = 0; i < n; ++i) if (DAG[i][j]) into[j]++;
	}
}



void toposort() {                    //拓扑排序
	cntx = 0;
	for (int i = 1; i <= n; ++i) {
		int j = 0;
		while (j < n && into[j] != 0) j++;
		Order[cntx++] = j + 1;
		into[j] = -inf;
		for (int k = 0; k < n; ++k) if (DAG[j][k]) into[k]--;
	}
}

bool dfs4(int u, int y, int z) {
	if (u == y) return true;
	vis2[u][0] = 1;
	for (int i = 0; i < n; i++)if (i != u && i != z && DAG[u][i] && !vis2[i][0]) {
		if (dfs4(i, y, z)) return true;
	}
	vis2[u][0] = 0;
	return false;
}

void dfs5(int d, int numm, int x, int y) {
	if (d >= numm) {
		int PAX1[10] = {}, PAX2[10] = {}, PAX3[10] = {}, cnt1 = 0, cnt2 = 0;
		for (int i = 0; i < numm; i++) {
			if (ess[i].second == 1) PAX1[cnt1++] = ess[i].first;
			else if (ess[i].second == 2) PAX2[cnt2++] = ess[i].first;
		}
		double tem1 = CH(LGObj, x, PAX1, cnt1);
		for (int i = 0; i < cnt2; i++)if (PAX2[i] != y) tem1 += CH(LGObj, PAX2[i], PAX3, 0);
		cnt1 = cnt2;
		for (int i = 0; i < numm; i++) {
			if (ess[i].second == 4) PAX1[cnt1++] = ess[i].first;
			else if (ess[i].second == 5) PAX2[cnt2++] = ess[i].first;
		}
		tem1 += CH(LGObj, y, PAX1, cnt1);
		for (int i = 0; i < cnt2; i++)if (PAX2[i] != x) tem1 += CH(LGObj, PAX2[i], PAX3, 0);
		if (tem1 > maxx) {
			maxx = tem1;
			for (int i = 0; i < numm; i++) {
				anss[i].first = ess[i].first;
				anss[i].second = ess[i].second;
			}
		}
		return;
	}
	if (ess[d].second == 3) {
		ess[d].second = 1;
		dfs5(d + 1, numm, x, y);
		ess[d].second = 2;
		dfs5(d + 1, numm, x, y);
		ess[d].second = 3;
	}
	else if (ess[d].second == 6) {
		ess[d].second = 4;
		dfs5(d + 1, numm, x, y);
		ess[d].second = 5;
		dfs5(d + 1, numm, x, y);
		ess[d].second = 6;
	}
	else dfs5(d + 1, numm, x, y);
}



void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {         //主函数
	int p, q, d, t, nC = 0, nF = 0, row, m, maxparent = 0, curnode;
	dlf *nodescore, lastscore, localmax, localscore = 0.0, *Temp;
	n = (int)mxGetScalar(mxGetField(prhs[0], 0, "VarNumber"));          //节点数
	LGObj.VarNumber = n;
	LGObj.CaseLength = (int)(mxGetScalar(mxGetField(prhs[0], 0, "CaseLength")));     //样本数量
	Temp = mxGetPr(mxGetField(prhs[0], 0, "VarRange"));           //每个节点可以取的值 例如1结点可以取值1 2 3那么VarRange[1][0]=1  VarRange[1][1]=2  VarRange[1][2]=3
	row = mxGetN(mxGetField(prhs[0], 0, "VarRange"));            //可以取值的数组行数

	for (p = 0; p < n; p++)
		for (q = 0; q < row; q++)
			LGObj.VarRange[p][q] = (int)Temp[p + q * n];
	for (int i = 0; i < n; i++) par[i] = i, rankk[i] = 0;

	Temp = mxGetPr(mxGetField(prhs[0], 0, "VarRangeLength"));      //每个节点可以取值的个数例如结点1可以取1 2 3那么值VarRangeLength[1]=3
	for (p = 0; p < n; p++) LGObj.VarRangeLength[p] = (int)Temp[p];
	Temp = mxGetPr(mxGetField(prhs[0], 0, "VarSample"));        //样本值

	for (p = 0; p < LGObj.CaseLength; p++)
		for (q = 0; q < n; q++) {
			LGObj.VarSample[p][q] = (int)Temp[p + q * LGObj.CaseLength];
		}
	Temp = mxGetPr(prhs[1]);
	for (p = 0; p < n; p++)
		for (q = 0; q < n; q++) {
			dag[p][q] = (int)Temp[p + q * n];
		}


	memset(DAG, 0, sizeof(DAG));
	memset(Order, 0, sizeof(Order));
	memset(vis3, 0, sizeof(vis3));
	memset(pnum1, 0, sizeof(pnum1));
	memset(pnum2, 0, sizeof(pnum2));
	memset(pnum3, 0, sizeof(pnum3));
	memset(vis, 0, sizeof(vis));
	memset(probilaty1, 0.0, sizeof(probilaty1));
	for (p = 0; p < n; p++) {                            //计算P(X=xi)(probilaty1)
		for (q = 0; q < LGObj.CaseLength; q++) {
			int k = LGObj.VarSample[q][p];
			probilaty1[p][k] = probilaty1[p][k] + 1.0;
		}
		for (int kk = 0; kk < LGObj.VarRangeLength[p]; kk++) {
			int k = LGObj.VarRange[p][kk];
			probilaty1[p][k] = (dlf)probilaty1[p][k] / LGObj.CaseLength;
			pnum1[p][k] = probilaty1[p][k];
		}
	}

	memset(probilaty2, 0.0, sizeof(probilaty2));
	memset(probilaty3, 0.0, sizeof(probilaty3));
	for (int p1 = 0; p1 < n; p1++)
		for (int p2 = 0; p2 < n; p2++) {
			for (int i = 0; i < LGObj.CaseLength; i++) {                          //计算P(X=xi,Y=yi)(probilaty2)和P(X=xi|Y=yi)(probilaty3)
				int k1 = LGObj.VarSample[i][p1], k2 = LGObj.VarSample[i][p2];
				probilaty2[p1][k1][p2][k2] = probilaty2[p1][k1][p2][k2] + 1.0;
			}
			for (int kk1 = 0; kk1 < LGObj.VarRangeLength[p1]; kk1++)
				for (int kk2 = 0; kk2 < LGObj.VarRangeLength[p2]; kk2++) {
					int k1 = LGObj.VarRange[p1][kk1], k2 = LGObj.VarRange[p2][kk2];
					probilaty2[p1][k1][p2][k2] = (dlf)probilaty2[p1][k1][p2][k2] / LGObj.CaseLength;
					pnum2[p1][k1][p2][k2] = probilaty2[p1][k1][p2][k2];
					if (probilaty1[p2][k2] >= eps) probilaty3[p1][k1][p2][k2] = probilaty2[p1][k1][p2][k2] / probilaty1[p2][k2];
				}
		}

	int cnt = 0;
	memset(MMI, 0.0, sizeof(MMI));
	for (p = 0; p < n; p++) {
		for (q = p + 1; q < n; q++) {                                      //计算 互信息(betentropy)
			dlf sum = 0.0;
			for (int kk1 = 0; kk1 < LGObj.VarRangeLength[p]; kk1++)
				for (int kk2 = 0; kk2 < LGObj.VarRangeLength[q]; kk2++) {
					int k1 = LGObj.VarRange[p][kk1], k2 = LGObj.VarRange[q][kk2];
					if (probilaty1[p][k1] >= eps && probilaty1[q][k2] >= eps && probilaty2[p][k1][q][k2] >= eps) sum += probilaty2[p][k1][q][k2] * log10(probilaty2[p][k1][q][k2] / (probilaty1[p][k1] * probilaty1[q][k2]));
				}
			betentropy[q][p] = betentropy[p][q] = sum;
			e[cnt++] = edge{ p,q,betentropy[p][q] };
		}
	}
	for (int i = 0; i < n; i++)
		for (int j = 0; j < n; j++)if (i != j) MMI[i] = max(MMI[i], betentropy[i][j]);  //计算每个结点最大的互信息


	for (int i = 0; i < n; i++) {
		P temm[maxs] = {};
		int cnt3 = 0;
		for (int j = 0; j < n; j++)if (i != j) temm[cnt3++] = P(betentropy[j][i], j);
		sort(temm, temm + cnt3);
		for (int j = cnt3 - 1; j >= 0 && j >= cnt3 - 3; j--) vis[i][temm[j].second] = 1;
	}

	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++)if (betentropy[i][j] >= alpha * MMI[i] && i != j) {     //满足条件(6)则连边
			unite(i, j);
			DAG[j][i] = DAG[i][j] = 1;
		}
	}

	sort(e, e + cnt);
	for (int i = 0; i < cnt; i++) {                  //将点对从大到小排序不在一个联通块的相连
		int u = e[i].u, v = e[i].v;
		if (found(u) == found(v)) continue;
		unite(u, v);
		DAG[u][v] = DAG[v][u] = 1;

	}



	memset(probilaty4, 0, sizeof(probilaty4));
	memset(probilaty5, 0, sizeof(probilaty5));               //计算P(X=xi,Y=yj,Z=zk)(probilaty4)和P(X=xi|Y=yj,Z=zk)(probilaty5)
	for (int i = 0; i < n; i++)
		for (int j = 0; j < n; j++)
			for (int k = 0; k < n; k++) {
				for (int h = 0; h < LGObj.CaseLength; h++) {
					int k1 = LGObj.VarSample[h][i], k2 = LGObj.VarSample[h][j], k3 = LGObj.VarSample[h][k];
					probilaty4[i][k1][j][k2][k][k3] = probilaty4[i][k1][j][k2][k][k3] + 1.0;
				}
				for (int kk1 = 0; kk1 < LGObj.VarRangeLength[i]; kk1++)
					for (int kk2 = 0; kk2 < LGObj.VarRangeLength[j]; kk2++)
						for (int kk3 = 0; kk3 < LGObj.VarRangeLength[k]; kk3++) {
							int k1 = LGObj.VarRange[i][kk1], k2 = LGObj.VarRange[j][kk2], k3 = LGObj.VarRange[k][kk3];
							probilaty4[i][k1][j][k2][k][k3] = (dlf)probilaty4[i][k1][j][k2][k][k3] / LGObj.CaseLength;
							pnum3[i][k1][j][k2][k][k3] = probilaty4[i][k1][j][k2][k][k3];
							if (probilaty2[j][k2][k][k3] >= eps) {
								probilaty5[i][k1][j][k2][k][k3] = probilaty4[i][k1][j][k2][k][k3] / probilaty2[j][k2][k][k3];
							}
						}
			}



	memset(e, 0, sizeof(e));
	cnt = 0;
	for (int i = 0; i < n; i++) {      //构造领接表
		P temm[maxs] = {};
		int cnt3 = 0; g[i].clear();
		for (int j = 0; j < n; j++)if (DAG[i][j] && i != j) temm[cnt3++] = P(betentropy[j][i], j);
		sort(temm, temm + cnt3);
		for (int j = cnt3 - 1; j >= 0; j--) {
			g[i].push_back(temm[j].second);
			//printf("%lf ", temm[j].first);
		}
		//printf("\n");
	}
	//printf("%lf %d %d\n", betentropy[2][7], vis[2][7], vis[7][2]);
	for (int i = 0; i < n; i++)
		for (int j = i + 1; j < n; j++)if (DAG[i][j]) e[cnt++] = edge{ i,j,betentropy[i][j] };
	sort(e, e + cnt);




	memset(viss, 0, sizeof(viss));             //找出三角环
	for (int i = 0; i < cnt; i++) {
		if (!DAG[e[i].u][e[i].v] || !DAG[e[i].v][e[i].u]) continue;
		for (int j = 0; j < g[e[i].v].size(); j++)if (DAG[e[i].u][g[e[i].v][j]] && DAG[g[e[i].v][j]][e[i].u] && g[e[i].v][j] != e[i].v) {
			int x = e[i].u, y = e[i].v, z = g[e[i].v][j];
			if (!DAG[y][z] || !DAG[z][y] || viss[x][y][z]) continue;
			double flag2 = 0;
			for (int kk1 = 0; kk1 < LGObj.VarRangeLength[x]; kk1++) {
				for (int kk2 = 0; kk2 < LGObj.VarRangeLength[y]; kk2++) {
					for (int kk3 = 0; kk3 < LGObj.VarRangeLength[z]; kk3++) {
						int k1 = LGObj.VarRange[x][kk1], k2 = LGObj.VarRange[y][kk2], k3 = LGObj.VarRange[z][kk3];
						if (fabs(pnum3[x][k1][y][k2][z][k3] * pnum1[z][k3]) >= eps && (pnum2[x][k1][z][k3] * pnum2[y][k2][z][k3]) >= eps)
							flag2 += pnum3[x][k1][y][k2][z][k3] * log(pnum3[x][k1][y][k2][z][k3] * pnum1[z][k3] / (pnum2[x][k1][z][k3] * pnum2[y][k2][z][k3]));
					}
				}
			}
			if (fabs(flag2) >= eps2) {
				int cnt2 = 3, PAX[maxs] = {};
				PAX[0] = x; PAX[1] = y; PAX[2] = z;
				for (int h = 0; h < n; h++)if (vis[h][x] && vis[h][y] && x != h && y != h && z != h) {//找出潜在的父节点
					PAX[cnt2++] = h;
					/*if (z == 5) printf("ac1 %d\n", h);
					if (x == 0 || y == 0) printf("ac2 %d\n", h)*/;
				}
				int flagg = 0;
				if (cnt2 > 3) {
					if (!calp(PAX, x, y, cnt2, 0)) {    //判断是否满足独立
						flagg = 1;
						break;
					}
				}
				if (!flagg && cnt2 > 3) {     //存在潜在的父节点并且满足条件独立
					DAG[x][y] = DAG[y][x] = 0;
					if (dag[x][y] || dag[y][x]) printf("wa1\n");
					for (int h = 0; h < cnt2; h++) { DAG[x][PAX[h]] = DAG[PAX[h]][x] = DAG[y][PAX[h]] = DAG[PAX[h]][y] = 1; }
					viss[x][y][z] = viss[x][z][y] = viss[y][x][z] = viss[y][z][x] = viss[z][x][y] = viss[z][y][x] = 1;
				}
			}
			else {          //满足循环则删除边x,y
				viss[x][y][z] = viss[x][z][y] = viss[y][x][z] = viss[y][z][x] = viss[z][x][y] = viss[z][y][x] = 1;
				DAG[x][y] = DAG[y][x] = 0;
				if (dag[x][y] || dag[y][x]) printf("wa1\n");
			}
		}
	}


	cnt = 0;
	for (int i = 0; i < n; i++) {      //更新领接表
		P temm[maxs] = {};
		int cnt3 = 0; g[i].clear();
		for (int j = 0; j < n; j++)if (DAG[i][j]) temm[cnt3++] = P(betentropy[i][j], j);
		sort(temm, temm + cnt3);
		for (int j = cnt3 - 1; j >= 0; j--) {
			g[i].push_back(temm[j].second);
		}
	}
	for (int i = 0; i < n; i++)
		for (int j = i + 1; j < n; j++)if (DAG[i][j]) e[cnt++] = edge{ i,j,betentropy[i][j] };
	sort(e, e + cnt);


	for (int i = 0; i < cnt; i++) {
		if (!DAG[e[i].u][e[i].v] || !DAG[e[i].v][e[i].u]) continue;
		for (int j = 0; j < g[e[i].v].size(); j++) {          //找出A-B-C不一定是三元环
			int x = e[i].u, z = e[i].v, y = g[e[i].v][j];
			if (!DAG[z][y] || !DAG[y][z] || y == x || dfs4(x, y, z) || DAG[x][y] || DAG[y][x]) continue;
			double flag2 = 0;
			for (int kk1 = 0; kk1 < LGObj.VarRangeLength[x]; kk1++) {
				for (int kk2 = 0; kk2 < LGObj.VarRangeLength[y]; kk2++) {
					for (int kk3 = 0; kk3 < LGObj.VarRangeLength[z]; kk3++) {
						int k1 = LGObj.VarRange[x][kk1], k2 = LGObj.VarRange[y][kk2], k3 = LGObj.VarRange[z][kk3];
						if ((pnum3[x][k1][y][k2][z][k3] * pnum1[z][k3]) >= eps && (pnum2[x][k1][z][k3] * pnum2[y][k2][z][k3]) >= eps)
							flag2 += 2 * pnum3[x][k1][y][k2][z][k3] * log(pnum3[x][k1][y][k2][z][k3] * pnum1[z][k3] / (pnum2[x][k1][z][k3] * pnum2[y][k2][z][k3]));
					}
				}
			}
			if (fabs(flag2) >= eps2) {
				int cnt2 = 3, PAX[maxs] = {};
				PAX[0] = x; PAX[1] = y; PAX[2] = z;
				//printf("%d %d %d\n", x, y, z);
				for (int h = 0; h < n; h++)if (vis[h][x] && vis[h][y] && x != h && y != h && z != h) {//添加父节点
					PAX[cnt2++] = h;
					//if (z == 5) printf("ac1 %d\n", h);
					//if (x == 0 || y == 0) printf("ac2 %d\n", h);
				}
				int flagg = 0;
				if (cnt2 > 3) {
					if (!calp(PAX, x, y, cnt2, 0)) {  //判断是否满足独立
						flagg = 1;
						break;
					}
				}
				else {                         //若潜在父节点为空则规定方向 A->C  B->C(3.3 phase4)
					printf("%lf\n", flag2);
					DAG[x][z] = DAG[y][z] = 1;
					DAG[z][x] = DAG[z][y] = 0;
					if (!dag[x][z] && dag[z][x]) printf("%d %d wa211\n", x, z);
					if (!dag[y][z] && dag[z][y]) printf("%d %d wa211\n", y, z);
					printf("ac211\n");
				}
				if (!flagg && cnt2 > 3) {
					DAG[x][y] = DAG[y][x] = 0;
					if (dag[x][y] || dag[y][x]) printf("wa111\n");
					for (int h = 0; h < cnt2; h++) { DAG[x][PAX[h]] = DAG[PAX[h]][x] = DAG[y][PAX[h]] = DAG[PAX[h]][y] = 1; }
				}
			}
			else {
				printf("ac1\n");
				if (dag[x][y] || dag[y][x]) printf("wa1\n");
				DAG[x][y] = DAG[y][x] = 0;
			}
		}
	}


	for (int i = 0; i < n; i++) {
		cnpa1 = 0; cnpa2 = 0;
		maxx = inf; cntx = 0; cnty = 0;      //cnty若为1则从i点出发存在四元环
		dfs2(i, i, -1, 0, 4);       //给四个结点组成的环评分
		if (cnty) {
			for (int k = 0; k < 3; k++) {
				int v1 = tem7[k], v2 = tem7[k + 1];
				if (tem6[k]) {
					DAG[v2][v1] = 0;
					DAG[v1][v2] = 1;
					if (dag[v2][v1] && !dag[v1][v2]) printf("%d %d wa3\n", v1, v2);
				}
				else {
					DAG[v1][v2] = 0;
					DAG[v2][v1] = 1;
					if (dag[v1][v2] && !dag[v2][v1]) printf("%d %d wa4\n", v1, v2);
				}
			}
			int v1 = tem7[3], v2 = i;
			if (tem6[3]) {
				DAG[v2][v1] = 0;
				DAG[v1][v2] = 1;
				if (dag[v2][v1] && !dag[v1][v2]) printf("%d %d wa5\n", v1, v2);
			}
			else {
				DAG[v1][v2] = 0;
				DAG[v2][v1] = 1;
				if (dag[v1][v2] && !dag[v2][v1]) printf("%d %d wa6\n", v1, v2);
			}
		}
	}



	for (int i = 0; i < n; i++) {       //更新领接表
		P temm[maxs] = {};
		int cnt3 = 0; g[i].clear();
		for (int j = 0; j < n; j++)if (DAG[i][j] && i != j) temm[cnt3++] = P(betentropy[i][j], j);
		sort(temm, temm + cnt3);
		for (int j = cnt3 - 1; j >= 0; j--) {
			g[i].push_back(temm[j].second);
		}
	}


	for (int i = 0; i < n; i++) {           //给五个结点组成的环评分
		cnpa1 = 0; cnpa2 = 0;
		maxx = inf; cntx = 0; cnty = 0;    //cnty若为1则从i点出发存在五元环
		dfs2(i, i, -1, 0, 5);
		if (cnty) {
			for (int k = 0; k < 4; k++) {
				int v1 = tem7[k], v2 = tem7[k + 1];
				if (tem6[k]) {
					DAG[v2][v1] = 0;
					DAG[v1][v2] = 1;
					if (dag[v2][v1] && !dag[v1][v2]) printf("%d %d wa7\n", v1, v2);
				}
				else {
					DAG[v1][v2] = 0;
					DAG[v2][v1] = 1;
					if (dag[v1][v2] && !dag[v2][v1]) printf("%d %d wa8\n", v1, v2);
				}
			}
			int v1 = tem7[4], v2 = i;
			if (tem6[4]) {
				DAG[v2][v1] = 0;
				DAG[v1][v2] = 1;
				if (dag[v2][v1] && !dag[v1][v2]) printf("%d %d wa9\n", v1, v2);
			}
			else {
				DAG[v1][v2] = 0;
				DAG[v2][v1] = 1;
				if (dag[v1][v2] && !dag[v2][v1]) printf("%d %d wa10\n", v1, v2);
			}
		}
	}

	for (int i = 0; i < n; i++) { g[i].clear(); fa[i].clear(); }  //faw记录结点i的所有父亲结点
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++)if (DAG[i][j]) {
			g[i].push_back(j); fa[j].push_back(i);
		}
	}


	P es[maxs*maxs] = {};
	int cntt = 0;
	for (int i = 0; i < n; i++) {
		int tem = 0;
		for (int j = 0; j < n; j++)if (i != j && DAG[j][i]) tem++;
		es[cntt++] = P(tem, i);
	}
	sort(es, es + cntt);


	for (int ii = cntt - 1; ii >= 0; ii--) {           //独立性检测(3.3 phase 6)找出x->y-z  x和y为有向  y到z是无向
		int i = es[ii].second;
		for (int j = 0; j < fa[i].size(); j++) {
			int v = fa[i][j];
			if (DAG[v][i] && !DAG[i][v]) {
				for (int k = 0; k < g[i].size(); k++)if (DAG[g[i][k]][i] && DAG[i][g[i][k]]) {
					int v2 = g[i][k];
					double flag2 = 0;
					if (DAG[v][v2] || DAG[v2][v]) continue;
					int x = v, y = i, z = v2;
					for (int kk1 = 0; kk1 < LGObj.VarRangeLength[x]; kk1++) {
						for (int kk2 = 0; kk2 < LGObj.VarRangeLength[y]; kk2++) {
							for (int kk3 = 0; kk3 < LGObj.VarRangeLength[z]; kk3++) {
								int k1 = LGObj.VarRange[x][kk1], k2 = LGObj.VarRange[y][kk2], k3 = LGObj.VarRange[z][kk3];
								if ((pnum3[x][k1][z][k3][y][k2] * pnum1[y][k2]) >= eps && (pnum2[x][k1][y][k2] * pnum2[z][k3][y][k2]) >= eps)
									flag2 += pnum3[x][k1][z][k3][y][k2] * log(pnum3[x][k1][z][k3][y][k2] * pnum1[y][k2] / (pnum2[x][k1][y][k2] * pnum2[z][k3][y][k2]));
							}
						}
					}
					if (fabs(flag2) < eps2) {
						DAG[i][v2] = DAG[v2][i] = 0;
						if (dag[i][v2] || dag[v2][i]) printf("wa20\n");
					}
					else {
						int cnt2 = 3, PAX[maxs] = {};
						PAX[0] = y; PAX[1] = z; PAX[2] = x;
						for (int h = 0; h < n; h++)if (vis[h][x] && vis[h][z] && x != h && y != h && z != h) {
							PAX[cnt2++] = h;
						}
						int flagg = 0;
						if (cnt2 > 3) {
							for (int i = 0; i < LGObj.VarRangeLength[y]; i++) {
								value[0] = LGObj.VarRange[y][i];
								if (!calp(PAX, x, z, cnt2, 1)) {
									flagg = 1;
									break;
								}
							}
						}
						else {
							DAG[v2][i] = 0; DAG[i][v2] = 1;
							printf("ac12\n");
							if (dag[v2][i]) printf("wa12\n");
							continue;
						}
						if (!flagg && cnt2 > 3) {
							DAG[x][z] = DAG[z][x] = 0;
							if (dag[x][z] || dag[z][x]) printf("wa19\n");
							for (int h = 0; h < cnt2; h++) { DAG[x][PAX[h]] = DAG[PAX[h]][x] = DAG[z][PAX[h]] = DAG[PAX[h]][z] = 1; }
						}
					}
				}
			}
		}
	}




	for (int i = 0; i < n; i++)      //(3.3 phase 7) 图切割
		for (int j = i + 1; j < n; j++)if (DAG[i][j] && DAG[j][i]) {
			int numm = 0;
			dlf temm1 = 0.0, temm2 = 0.0;
			for (int k = 0; k < n && numm < 4; k++)if ((DAG[k][i] || DAG[i][k]) && k != j) {
				if (DAG[k][i] && DAG[i][k]) ess[numm++] = P(k, 3);
				else if (DAG[k][i]) ess[numm++] = P(k, 1);
				else ess[numm++] = P(k, 2);

			}
			ess[numm++] = P(j, 3);
			for (int k = 0; k < n && numm < 5; k++)if ((DAG[k][j] || DAG[j][k]) && k != i) {
				if (DAG[k][j] && DAG[j][k]) ess[numm++] = P(k, 6);
				else if (DAG[k][j]) ess[numm++] = P(k, 4);
				else ess[numm++] = P(k, 5);
			}
			maxx = inf;
			dfs5(0, numm, i, j);
			for (int k = 0; k < numm; k++)if (ess[k].second == 3 || ess[k].second == 6) {
				int u = anss[k].first;
				if (anss[k].second == 1) {
					DAG[i][u] = 1; DAG[u][i] = 0;
					if (!dag[i][u] && dag[u][i]) printf("%d %d wa15\n", u, i);
					printf("ac15\n");
				}
				else if (anss[k].second == 2) {
					DAG[u][i] = 1; DAG[i][u] = 0;
					if (!dag[u][i] && dag[i][u]) printf("%d %d wa16\n", u, i);
					printf("ac16\n");
				}
				else if (anss[k].second == 4) {
					DAG[j][u] = 1; DAG[u][j] = 0;
					if (!dag[j][u] && dag[u][j]) printf("%d %d wa17\n", u, j);
					printf("ac17\n");
				}
				else {
					DAG[u][j] = 0; DAG[j][u] = 1;
					if (dag[u][j] && !dag[j][u]) printf("%d %d wa18\n", u, j);
					printf("ac18\n");
				}
			}

			for (int k = 0; k < n; k++) fa[k].clear();
			for (int k = 0; k < n; k++)
				for (int kk = 0; kk < n; kk++)if (DAG[kk][k])  fa[k].push_back(kk);

		}

	plhs[0] = mxCreateDoubleMatrix(n, n, mxREAL);
	Temp = mxGetPr(plhs[0]);
	for (p = 0; p < n; p++)
		for (q = 0; q < n; q++)
			Temp[p + q * n] = (dlf)DAG[p][q];


	get_into_degree();    //拓扑排序
	toposort();           //拓扑排序按字典序排
	plhs[1] = mxCreateDoubleMatrix(1, n, mxREAL);
	Temp = mxGetPr(plhs[1]);
	for (int i = 0; i < cntx; i++) Temp[i] = (dlf)Order[i];
}






