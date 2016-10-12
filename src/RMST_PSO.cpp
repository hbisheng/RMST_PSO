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

/*********************************** �㷨����  	*****************************************/
const int N_particle = 50;  	//���Ӹ���
const int N_MAXITER	 = 1000;	//��������
const double c1 = 1.0 ;			//ѧϰ����1 
const double c2 = 1.0 ;			//ѧϰ����2 

 
/*********************************** ���߷�Χ��Ŀ��㣬Hanan�� **************************/
const int WEIGHT = 10000;
const int HEIGHT = 10000;

const int N_maxpoint = 999; 	//���������������� ��Ҫ���� N_origin + N_origin - 2 
const int N_maxhanan = 9999; 	//Hanan�������������Ҫ���� N_origin * N_origin - 1 

int N_origin 	= 0; 	//�������
int N_hanan 	= 0;	//Hanan����� 
pair<int, int> pos_origin[N_maxpoint];	//��������� 
pair<int, int> pos_hanan[N_maxhanan];	//Hanan������ 
map<pair<int, int>, int> map_point; 	//�㵽intֵ�����䣬����ȥ�� 

/*********************************** PSO������Ϣ  *****************************************/
int N_steiner_particle[N_particle] = {0};  //�����ӵ�steiner����� 
pair<double,double> v_particle[N_particle][N_maxpoint] ;	//��ǰ�ٶ� 
pair<double,double> pos_particle[N_particle][N_maxpoint]; 	//��ǰλ�� 
pair<int, int>		edge_particle[N_particle][N_maxpoint];	//��С�������߱� 

//����������Ϣ 
int value_pbest[N_particle] 	= {0};		//���Ž��ֵ 
int N_steiner_pbest[N_particle] = {0};		//���Ž��steiner����� 
pair<double,double> pos_pbest[N_particle][N_maxpoint];	//���Ž������λ�� 
pair<int, int> edge_pbest[N_particle][N_maxpoint];		//���Ž�ı߱� 

//ȫ��������Ϣ 
int value_gbest = 0x7FFFFFFF;	//���Ž��ֵ
int N_steiner_gbest = 0;	//���Ž��steiner����� 
pair<double,double> pos_gbest[N_maxpoint]; 	//���Ž������λ�� 
pair<int, int> edge_gbest[N_maxpoint];		//���Ž�ı߱� 
/*********************************** ��С��������غ���  *****************************************/

int DIS(pair<int,int> a, pair<int,int> b) //������֮��ľ��� 
{
	return abs(a.first - b.first) + abs(a.second - b.second);
}
pair<int,int> get_hanan(pair<double,double> s) //��ȡ��������λ�������hanan�� 
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

/*********************************** �����С������  *****************************************/
int mst_mark[N_maxpoint];
pair<int, int> mst_point[N_maxpoint];
int get_mst(int particle) //�����С���������� 
{	
	//�ȼ���ԭʼ�� 
	int mst_num = 0;
	for(int i = 0; i != N_origin; i++)
		mst_point[mst_num++] = pos_origin[i];
	
	//�ټ���steiner�� 
	int N_candidate = N_steiner_particle[particle];
	map_point.clear();
	map<pair<int,int>,int>::iterator iter;
	for(int i = 0; i < N_candidate; i++)
	{
		pair<int, int> candidate = get_hanan(pos_particle[particle][i]);//��ȡ������iά���������hanan��
		iter = map_point.find(candidate);
		if(iter == map_point.end())
		{
			mst_point[mst_num++] = candidate;
			map_point[candidate] = 1;//����� 
		}
	}
	
	//PRIM��С�������㷨 
	memset(mst_mark, 0, sizeof(mst_mark));
	int u_num = 0,		//��ǰ�������ĵ��� 
		edge_num = 0,	//��ǰ�ߵ����� 
		sum_len = 0;	//��ǰ�߳�
			 
	mst_mark[0] = 1; // ȡ��һ������Ϊ��� 
	u_num++;
	while(u_num < mst_num) //�����е�δ���뵱ǰ������ʱ 
	{
		int min_idx = 0, min_len = 0x7FFFFFFF;
		for(int i = 0; i < mst_num; i++)
			for(int j = i + 1; j < mst_num; j++)
				if(mst_mark[i] + mst_mark[j] == 1)//һ�����У�һ���������� 
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

/*********************************** ���¾ֲ���ȫ������ֵ  *****************************************/
void refresh_pbest(int particle, int value)
{
	N_steiner_pbest[particle] = N_steiner_particle[particle];
	value_pbest[particle] = value;
	for(int i = 0; i < N_steiner_particle[particle]; i++)
		pos_pbest[particle][i] = pos_particle[particle][i];
		
	for(int i = 0; i < N_steiner_particle[particle] + N_origin - 1; i++) //�� = ���� - 1 
		edge_pbest[particle][i] = edge_particle[particle][i];				
}
void refresh_gbest(int particle, int value) 
{
	N_steiner_gbest = N_steiner_particle[particle];
	value_gbest = value;
	for(int i = 0; i < N_steiner_particle[particle]; i++)
		pos_gbest[i] = pos_particle[particle][i];
	for(int i = 0; i < N_steiner_particle[particle] + N_origin  - 1; i++) //�� = ���� - 1 
		edge_gbest[i] = edge_particle[particle][i];
}	

/*********************************** ���������Ľ�PSO�������RMST  *********************************/
int main()
{
	//Step0���������� 
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
		map_point[pos_origin[i]] = 1; //�е��λ�����ϱ�� 
	}
	cerr << "��С���������ȣ� " << get_mst(0) << endl;  
	
	
	
	//Step1������Hanan��
	map<pair<int,int>,int>::iterator iter; 
	for(int p_one = 0; p_one < N_origin; p_one++)
		for(int p_two = p_one + 1; p_two < N_origin; p_two++)
			if(pos_origin[p_one].first != pos_origin[p_two].first && pos_origin[p_one].second != pos_origin[p_two].second) //������ͬһ�л����� 
			{
				pair<int, int> p1 = make_pair(pos_origin[p_one].first,pos_origin[p_two].second); //һ�����ܵ���Hanan�� 
				iter = map_point.find(p1);
				if(iter == map_point.end())//�����δ������ͼ�� 
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
		cout << "ԭͼ����С��������Ϊ����RMST" << endl;
		return 0;
	}
	cout << "--Step2" << endl;
	//Step2����ʼ������λ�á��ٶȵ� 
	srand(time(NULL));
	int N_steiner_limit = min(N_origin - 2, N_hanan);	//steiner�������������� 
	//��ʼ��һ���������
	int random_order[N_maxpoint] = {0};
	for(int i = 0; i != N_steiner_limit; i++)
		random_order[i] = i;
	for(int particle = 0; particle != N_particle; particle++)//��ÿһ�����ӳ�ʼ��
	{
		//�ٶ�v_particle�Ѿ���ʼ��Ϊ0
		//����ά����һ�����hanan����Ϊ��ʼλ�ã�
		for(int n = 0; n != N_steiner_limit; n++)
			swap(random_order[rand()%N_steiner_limit], random_order[rand()%N_steiner_limit]);	
		for(int k = 0; k != N_steiner_limit; k++)
		 	pos_particle[particle][k] = pos_hanan[random_order[k]];
		//steiner����[0,n_steiner_limit]��� 
		N_steiner_particle[particle] = rand() % (N_steiner_limit+1);
		
		//��ʼ����������ֵ 
		int plen = get_mst(particle);
		
		refresh_pbest(particle, plen);
		//����ȫ������ֵ 
		if(plen < value_gbest)
			refresh_gbest(particle, plen);
	}
	cout << "--Step3  " << value_gbest << endl;
	//Step3: ��ʼ����
	int n_iter = 0;
	int best_last = value_gbest;
	while(n_iter++ < N_MAXITER)
	{
		//���¸������ӵ��ٶȺ�λ�� 
		double w = 0.9 - n_iter * 0.5 / N_MAXITER;//����ϵ�� 
		for(int particle = 0; particle != N_particle; particle++)
		{ 
			for(int dim = 0; dim != N_steiner_particle[particle]; dim++)
			{
				//r1, r2��(0,1)������С�� 
				double 	r1 = (rand() % 30000) / 30000.0,
				 		r2 = (rand() % 30000) / 30000.0;
				//��������dimά���ٶ� 
				v_particle[particle][dim].first = w * v_particle[particle][dim].first 
										+ c1 * r1 * (pos_pbest[particle][dim].first - pos_particle[particle][dim].first)
										+ c2 * r2 * (pos_gbest[dim].first - pos_particle[particle][dim].first);  
				v_particle[particle][dim].second = w * v_particle[particle][dim].second 
										+ c1 * r1 * (pos_pbest[particle][dim].second - pos_particle[particle][dim].second)
										+ c2 * r2 * (pos_gbest[dim].second - pos_particle[particle][dim].second);  
				//��������dimά��λ��
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
			//����steiner����� 
			
			int new_steiner = max(max(N_steiner_particle[particle], N_steiner_pbest[particle]), N_steiner_gbest); 
			//��������ˣ�Ӧ�ö��¿���ά�ȵ�pbest�����޸ģ���Ȼ�ᵽ�����ܡ� 
			if(new_steiner != N_steiner_particle[particle])
			{
				for(int n_steiner = N_steiner_particle[particle]; n_steiner != new_steiner; n_steiner++)//�����ӵ�ά��Ҫ��ʼ��pbest 
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
	
	
	cout << "------------------------- ��������λ�� ------------------------------" << endl;
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
	//���ȫ������ֵ
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
