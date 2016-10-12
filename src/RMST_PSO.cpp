#include<iostream>
#include<cmath>
#include<cstdlib>
#include<cstdio>
#include<ctime>
#include<utility>
#include<memory.h>
#include<map>
#include<float.h>

#define max(x,y) ((x)>(y)?(x):(y))
#define min(x,y) ((x)<(y)?(x):(y))
using namespace std;

/*********************************** 算法参数  	*****************************************/
const int N_particle = 50;  	//粒子个数
const int N_MAXITER	 = 1000;	//迭代次数
const double c1 = 1.0 ;			//学习参数1 
const double c2 = 1.0 ;			//学习参数2 

 
/*********************************** 布线范围，目标点，Hanan点 **************************/
const int WEIGHT = 10000;
const int HEIGHT = 10000;

const int N_maxpoint = 999; 	//生成树的最大点数， 需要大于 N_origin + N_origin - 2 
const int N_maxhanan = 9999; 	//Hanan点的最大点数，需要大于 N_origin * N_origin - 1 

int N_origin 	= 0; 	//输入点数
int N_hanan 	= 0;	//Hanan点个数 
pair<int, int> pos_origin[N_maxpoint];	//输入点坐标 
pair<int, int> pos_hanan[N_maxhanan];	//Hanan点坐标 
map<pair<int, int>, int> map_point; 	//点到int值的隐射，用于去重 

/*********************************** PSO粒子信息  *****************************************/
int N_steiner_particle[N_particle] = {0};  //各粒子的steiner点个数 
pair<double,double> v_particle[N_particle][N_maxpoint] ;	//当前速度 
pair<double,double> pos_particle[N_particle][N_maxpoint]; 	//当前位置 
pair<int, int>		edge_particle[N_particle][N_maxpoint];	//最小生成树边表 

//个体最优信息 
int value_pbest[N_particle] 	= {0};		//最优解的值 
int N_steiner_pbest[N_particle] = {0};		//最优解的steiner点个数 
pair<double,double> pos_pbest[N_particle][N_maxpoint];	//最优解的粒子位置 
pair<int, int> edge_pbest[N_particle][N_maxpoint];		//最优解的边表 

//全局最优信息 
int value_gbest = 0x7FFFFFFF;	//最优解的值
int N_steiner_gbest = 0;	//最优解的steiner点个数 
pair<double,double> pos_gbest[N_maxpoint]; 	//最优解的粒子位置 
pair<int, int> edge_gbest[N_maxpoint];		//最优解的边表 
/*********************************** 最小生成树相关函数  *****************************************/

int DIS(pair<int,int> a, pair<int,int> b) //两个点之间的距离 
{
	return abs(a.first - b.first) + abs(a.second - b.second);
}
pair<int,int> get_hanan(pair<double,double> s) //获取距离粒子位置最近的hanan点 
{
	int min_idx = -1; 
	double min_len = DBL_MAX;
	for(int i = 0; i != N_hanan; i++)
	{
		double plen = abs(s.first - pos_hanan[i].first) + abs(s.second - pos_hanan[i].second);
		if(plen < min_len)
		{
			min_len = plen;
			min_idx = i;
		}
	}
	return pos_hanan[min_idx];
}

/*********************************** 求解最小生成树  *****************************************/
int mst_mark[N_maxpoint];
pair<int, int> mst_point[N_maxpoint];
int get_mst(int particle) //求解最小生成树长度 
{	
	//先加入原始点 
	int mst_num = 0;
	for(int i = 0; i != N_origin; i++)
		mst_point[mst_num++] = pos_origin[i];
	
	//再加入steiner点 
	int N_candidate = N_steiner_particle[particle];
	map_point.clear();
	map<pair<int,int>,int>::iterator iter;
	for(int i = 0; i < N_candidate; i++)
	{
		pair<int, int> candidate = get_hanan(pos_particle[particle][i]);//获取粒子在i维度上最近的hanan点
		iter = map_point.find(candidate);
		if(iter == map_point.end())
		{
			mst_point[mst_num++] = candidate;
			map_point[candidate] = 1;//做标记 
		}
	}
	
	//PRIM最小生成树算法 
	memset(mst_mark, 0, sizeof(mst_mark));
	int u_num = 0,		//当前生成树的点数 
		edge_num = 0,	//当前边的条数 
		sum_len = 0;	//当前边长
			 
	mst_mark[0] = 1; // 取第一个点作为起点 
	u_num++;
	while(u_num < mst_num) //当还有点未加入当前生成树时 
	{
		int min_idx = 0, min_len = 0x7FFFFFFF;
		for(int i = 0; i < mst_num; i++)
			for(int j = i + 1; j < mst_num; j++)
				if(mst_mark[i] + mst_mark[j] == 1)//一个树中，一个树外的情况 
					if( DIS(mst_point[i], mst_point[j]) < min_len)
					{
						min_len = DIS(mst_point[i], mst_point[j]);
						if(mst_mark[i]) 
							min_idx = j;
						else
							min_idx = i;
						edge_particle[particle][edge_num] = make_pair(i,j); 
					}
		mst_mark[min_idx] = 1;
		u_num++;
		edge_num ++;
		sum_len += min_len;
	}
	//cout << "U E " << u_num << ' ' << edge_num << endl;
	return sum_len;
}

/*********************************** 更新局部、全局最优值  *****************************************/
void refresh_pbest(int particle, int value)
{
	N_steiner_pbest[particle] = N_steiner_particle[particle];
	value_pbest[particle] = value;
	for(int i = 0; i < N_steiner_particle[particle]; i++)
		pos_pbest[particle][i] = pos_particle[particle][i];
		
	for(int i = 0; i < N_steiner_particle[particle] + N_origin - 1; i++) //边 = 点数 - 1 
		edge_pbest[particle][i] = edge_particle[particle][i];				
}
void refresh_gbest(int particle, int value) 
{
	N_steiner_gbest = N_steiner_particle[particle];
	value_gbest = value;
	for(int i = 0; i < N_steiner_particle[particle]; i++)
		pos_gbest[i] = pos_particle[particle][i];
	for(int i = 0; i < N_steiner_particle[particle] + N_origin  - 1; i++) //边 = 点数 - 1 
		edge_gbest[i] = edge_particle[particle][i];
}	

/*********************************** 主函数：改进PSO方法求解RMST  *********************************/
int main()
{
	//Step0：读入数据 
	freopen("input.txt", "r", stdin);
	freopen("ans.txt", "w", stdout);
	cin	>> N_origin;
	for(int i = 0; i != N_origin; i++)
	{
		cin >> pos_origin[i].first >> pos_origin[i].second;		
		if(pos_origin[i].first > WEIGHT || pos_origin[i].second > HEIGHT || pos_origin[i].first < 0 || pos_origin[i].second < 0)
		{
			cout << "Index out of bound!" << endl;
			return -1; 
		}
		map_point[pos_origin[i]] = 1; //有点的位置做上标记 
	}
	cerr << "最小生成树长度： " << get_mst(0) << endl;  
	
	
	
	//Step1：生成Hanan点
	map<pair<int,int>,int>::iterator iter; 
	for(int p_one = 0; p_one < N_origin; p_one++)
		for(int p_two = p_one + 1; p_two < N_origin; p_two++)
			if(pos_origin[p_one].first != pos_origin[p_two].first && pos_origin[p_one].second != pos_origin[p_two].second) //不能在同一行或列上 
			{
				pair<int, int> p1 = make_pair(pos_origin[p_one].first,pos_origin[p_two].second); //一个可能的新Hanan点 
				iter = map_point.find(p1);
				if(iter == map_point.end())//如果还未存在于图上 
				{
					pos_hanan[N_hanan++] = p1;
					map_point[p1] = 1;
				}

				pair<int, int> p2 = make_pair(pos_origin[p_two].first,pos_origin[p_one].second);
				iter = map_point.find(p2);
				if(iter == map_point.end())
				{
					pos_hanan[N_hanan++] = p2;
					map_point[p2] = 1;
				}
			}
	cout << "N_hanan: " << N_hanan << endl;
	if(!N_hanan)
	{
		cout << "原图的最小生成树即为所求RMST" << endl;
		return 0;
	}
	cout << "--Step2" << endl;
	//Step2：初始化粒子位置、速度等 
	srand(time(NULL));
	int N_steiner_limit = min(N_origin - 2, N_hanan);	//steiner点的最大数量限制 
	//初始化一个随机排列
	int random_order[N_maxpoint] = {0};
	for(int i = 0; i != N_steiner_limit; i++)
		random_order[i] = i;
	for(int particle = 0; particle != N_particle; particle++)//对每一个粒子初始化
	{
		//速度v_particle已经初始化为0
		//各个维度以一个随机hanan点作为初始位置，
		for(int n = 0; n != N_steiner_limit; n++)
			swap(random_order[rand()%N_steiner_limit], random_order[rand()%N_steiner_limit]);	
		for(int k = 0; k != N_steiner_limit; k++)
		 	pos_particle[particle][k] = pos_hanan[random_order[k]];
		//steiner点在[0,n_steiner_limit]随机 
		N_steiner_particle[particle] = rand() % (N_steiner_limit+1);
		
		//初始化个体最优值 
		int plen = get_mst(particle);
		
		refresh_pbest(particle, plen);
		//更新全局最优值 
		if(plen < value_gbest)
			refresh_gbest(particle, plen);
	}
	cout << "--Step3  " << value_gbest << endl;
	//Step3: 开始迭代
	int n_iter = 0;
	int best_last = value_gbest;
	while(n_iter++ < N_MAXITER)
	{
		//更新各个粒子的速度和位置 
		double w = 0.9 - n_iter * 0.5 / N_MAXITER;//惯性系数 
		for(int particle = 0; particle != N_particle; particle++)
		{ 
			for(int dim = 0; dim != N_steiner_particle[particle]; dim++)
			{
				//r1, r2是(0,1)间的随机小数 
				double 	r1 = (rand() % 30000) / 30000.0,
				 		r2 = (rand() % 30000) / 30000.0;
				//更新粒子dim维的速度 
				v_particle[particle][dim].first = w * v_particle[particle][dim].first 
										+ c1 * r1 * (pos_pbest[particle][dim].first - pos_particle[particle][dim].first)
										+ c2 * r2 * (pos_gbest[dim].first - pos_particle[particle][dim].first);  
				v_particle[particle][dim].second = w * v_particle[particle][dim].second 
										+ c1 * r1 * (pos_pbest[particle][dim].second - pos_particle[particle][dim].second)
										+ c2 * r2 * (pos_gbest[dim].second - pos_particle[particle][dim].second);  
				//更新粒子dim维的位置
				pos_particle[particle][dim].first 	+= v_particle[particle][dim].first; 
				pos_particle[particle][dim].second 	+= v_particle[particle][dim].second; 
			}	
			int plen = get_mst(particle);
			if(plen < value_pbest[particle])
			{
				cout << "PPP_refresh, n_iter: " << n_iter << " particle: " << particle << endl; 	
				refresh_pbest(particle,plen);
			}
			if(plen < value_gbest)
			{
				cout << "GGGGGGG_refresh, n_iter: " << n_iter << " particle: " << particle << endl; 	
				refresh_gbest(particle,plen);
			}
			//更新steiner点个数 
			
			int new_steiner = max(max(N_steiner_particle[particle], N_steiner_pbest[particle]), N_steiner_gbest); 
			//如果更新了，应该对新开的维度的pbest进行修改，不然会到处乱跑。 
			if(new_steiner != N_steiner_particle[particle])
			{
				for(int n_steiner = N_steiner_particle[particle]; n_steiner != new_steiner; n_steiner++)//新增加的维度要初始化pbest 
					pos_pbest[particle][n_steiner] = pos_particle[particle][n_steiner];
				N_steiner_particle[particle] = new_steiner;
			}
		} 
		
		int p = 0;
		cout << "Particle: " << p << endl; 
		for(int dim = 0; dim != N_steiner_particle[p]; dim++)
		{
			pair<int, int> tmp = get_hanan(pos_particle[p][dim]);
			cout << "----- DIM: " << dim << ' ' 
				 << int(pos_particle[p][dim].first) << ' ' << int(pos_particle[p][dim].second) 
				 << " TRANSLATED: " << tmp.first << ' ' << tmp.second
				 << " LAST_VELOCITY: " << int(v_particle[p][dim].first) << ' ' << int(v_particle[p][dim].second) 
				 << " PBEST: " << int(pos_pbest[p][dim].first) << ' ' << int(pos_pbest[p][dim].second)  
				 << " GBEST: " << int(pos_gbest[dim].first) << ' ' << int(pos_gbest[dim].second) << endl;
		}
		
		/* 
		if(value_gbest != best_last)
		{
			cout << n_iter << ' ' << best_last << " -> " << value_gbest << endl; 
			best_last = value_gbest;
		}
		*/
	} 
	
	
	cout << "------------------------- 最后的粒子位置 ------------------------------" << endl;
	for(int particle = 0; particle != N_particle; particle++)
	{
		cout << "Particle: " << particle << endl; 
		for(int dim = 0; dim != N_steiner_particle[particle]; dim++)
		{
			pair<int, int> tmp = get_hanan(pos_particle[particle][dim]);
			cout << "----- DIM: " << dim << ' ' 
				 << pos_particle[particle][dim].first << ' ' << pos_particle[particle][dim].second 
				 << " TRANSLATED: " << tmp.first << ' ' << tmp.second << endl; 
		}
	}
	
	
	cout << "----------------------------------------------------------------------" << endl; 
	//输出全局最优值
	cout << "Steiner Points:" << endl;
	for(int i = 0; i != N_steiner_gbest; i++)
	{ 
		pair<int, int> p = get_hanan(pos_gbest[i]);
		cout << p.first << ' ' << p.second << endl;
	} 
	cout << "Edges:" << endl;
	for(int i = 0; i != N_steiner_gbest + N_origin - 1; i++)
		cout << edge_gbest[i].first << ' ' << edge_gbest[i].second << endl;
	cerr << "Value:" << value_gbest << endl;
	
	cerr << "Time:" << double(clock()) / CLOCKS_PER_SEC << endl;
	return 0;
} 
