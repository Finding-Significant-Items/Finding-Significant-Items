/*
统计输出的次数
*/
#include <iostream>
#include <memory.h>
#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <fstream>
#include <functional>
#include <ctime>
#include <map>
#include <cstring>
#include <string>
#include <algorithm>
#include <set>


using namespace std;

//----------------------------------the parameters of STBF

double ARE = 0;
double hash_A = (sqrt(5) - 1) / 2;
double hash_B = (sqrt(7) - 2) / 2;
double hash_C = sqrt(3) / 2;

// mem = 4*a*m
// change k , a , m
int a,b,m;  // a means the number of peirods,  b means the number of items in each peiord, m means the numer of cells in each STBF
#define w1 9600000	//the number of counters
#define hash_print_max 0xFFFF	//the maximum value of hash finger
#define hash_print_1_max 0xFF	//the another hash finger
#define int_max 0xFFFFFFFF	
#define g 104				//Let threshold g=bit_of_ID / bit_of_Raptor
#define bit_of_ID 32	//the length of an item ID
#define str_max_len 4
#define STR_MAX_LEN_INPUT 4



#define th1 3	
#define hash_number	3



//4 bytes
class STBF {
public:
	bool raptor;
	unsigned short hash_print;
	unsigned char hash_print_1;
	bool flag;
	STBF() {
		raptor = 0;
		hash_print = 0;
		flag = 0;
		hash_print_1 = 0;
	};
};


class ID {
public:
	char x[str_max_len];
	bool if_equal(ID n) {
		bool flag = 1;
		for (int i = 0; i < str_max_len;i++) {
			if (x[i] != n.x[i]) {
				flag = 0;
				break;
			}
		}
		return flag;
	}
	ID() {
		memset(x, 0, sizeof(x));
	}
};

class ID_input {
public:
	char x[STR_MAX_LEN_INPUT];
	ID_input() {
		memset(x, 0, sizeof(x));
	}
};

bool operator < (ID an, ID bn) {
	for (int i = 0;i < str_max_len;i++) {
		if (bn.x[i] < an.x[i]) {
			return true;
		}
		else if (bn.x[i] > an.x[i]) {
			return false;
		}
	}
	return false;
}

ID_input * data_of_id;
int * Time;
unsigned char L1[w1] = { 0 };
bool hash_print_memory[hash_print_max] = { 0 };	//to record whether a hash finger has tried to be decoded
STBF *s; // STBF for hot
multimap<ID, int> id_map;
map<ID, int> all_id_map;
set<ID> id_set;
unsigned int murmur3_32(const char* key, size_t len, unsigned int seed) {
	unsigned int h = seed;
	if (len > 3) {
		const unsigned int* key_x4 = (const unsigned int*)key;
		size_t i = len >> 2;
		do {
			unsigned int k = *key_x4++;
			k *= 0xcc9e2d51;
			k = (k << 15) | (k >> 17);
			k *= 0x1b873593;
			h ^= k;
			h = (h << 13) | (h >> 19);
			h = (h * 5) + 0xe6546b64;
		} while (--i);
		key = (const char*)key_x4;
	}
	if (len & 3) {
		size_t i = len & 3;
		unsigned int k = 0;
		key = &key[i - 1];
		do {
			k <<= 8;
			k |= *key--;
		} while (--i);
		k *= 0xcc9e2d51;
		k = (k << 15) | (k >> 17);
		k *= 0x1b873593;
		h ^= k;
	}
	h ^= len;
	h ^= h >> 16;
	h *= 0x85ebca6b;
	h ^= h >> 13;
	h *= 0xc2b2ae35;
	h ^= h >> 16;
	return h;
}

int hash_1(int hash_range, unsigned int ID) {
	double x = ID * hash_A;
	double y = x - (unsigned int)x;
	int z = (int)(hash_range * y);
	return z;
};

int hash_2(int hash_range, unsigned int ID) {
	double x = ID * hash_B;
	double y = x - (unsigned int)x;
	int z = (int)(hash_range * y);
	return z;
}

int hash_3(int hash_range, unsigned int ID) {
	double x = ID * hash_C;
	double y = x - (unsigned int)x;
	int z = (int)(hash_range * y);
	return z;
}

unsigned short hash_o(unsigned int ID) {
	double x = ID * hash_A;
	double y = x - (unsigned int)x;
	unsigned short z = (unsigned short)(hash_print_max * y);
	return z;
}

unsigned char hash_o_1(unsigned int ID) {
	double x = ID * hash_B;
	double y = x - (unsigned int)x;
	unsigned char z = (unsigned char)(hash_print_1_max * y);
	return z;
}

//1 for hot and 0 for cold
bool cold_filter(unsigned int ID) {
	unsigned char V1 = th1;
	int h_L1[hash_number];
	h_L1[0] = hash_1(w1, ID);
	h_L1[1] = hash_2(w1, ID);
	h_L1[2] = hash_3(w1, ID);
	for (int j = 0; j < hash_number; j++) {
		if (L1[h_L1[j]] < V1) {
			V1 = L1[h_L1[j]];
		}
	}
	if (V1 < th1) {
		for (int j = 0; j < hash_number; j++) {
			if (L1[h_L1[j]] == V1) {
				L1[h_L1[j]] = L1[h_L1[j]] + 1;
			}
		}
		return 0;
	}
	else {
		return 1;
	}
}

bool encode(ID myid, int seed) {
	srand(seed);
	bool coeff_a[bit_of_ID];
	for (int i = 0; i < bit_of_ID;i++) {
		coeff_a[i] = rand() % 2;
	}



	bool ID_binary[bit_of_ID];
	for (int i = 0; i < str_max_len; i++) {
		for (int j = 0; j < 8;j++) {
			ID_binary[i * 8 + j] = bool(myid.x[i] & (1 << (8 - 1 - j)));
		}
	}
	bool Raptor = 0;
	for (int i = 0; i < bit_of_ID; i++) {
		Raptor = Raptor ^ (coeff_a[i] & ID_binary[i]);
	}

	return Raptor;
}

//Gauss
ID gauss(int _cnt, bool long_coeff[][bit_of_ID], bool *long_raptor) {


	for (int i = 0; i < bit_of_ID; i++) {
		if (long_coeff[i][i] == 0) {
			for (int j = i + 1; j < _cnt; j++) {
				if (long_coeff[j][i] == 1) {
					for (int k = 0; k < bit_of_ID; k++) {
						long_coeff[i][k] = long_coeff[i][k] ^ long_coeff[j][k];
					}
					long_raptor[i] = long_raptor[i] ^ long_raptor[j];
					break;
				}
			}
		}

		for (int j = i + 1; j < _cnt; j++) {
			if (long_coeff[j][i] == 1) {
				for (int k = 0; k < bit_of_ID; k++) {
					long_coeff[j][k] = long_coeff[j][k] ^ long_coeff[i][k];
				}
				long_raptor[j] = long_raptor[j] ^ long_raptor[i];
			}

		}

	}



	for (int i = bit_of_ID - 1; i >= 0; i--) {
		for (int j = i - 1; j >= 0; j--) {
			if (long_coeff[j][i] == 1) {
				long_coeff[j][i] = long_coeff[j][i] ^ long_coeff[i][i];
				long_raptor[j] = long_raptor[j] ^ long_raptor[i];
			}
		}
	}



	ID myid;


	for (int i = 0; i < str_max_len; i++) {
		for (int j = 0; j < 8;j++) {
			if (long_raptor[i * 8 + j]) {
				myid.x[i] = myid.x[i] | (1 << (8 - 1 - j));
			}
		}
	}

	return myid;
}



//construct coefficients
void parameter(int &_cnt, int *cnt, bool *raptor, bool  long_coeff[][bit_of_ID], bool * long_raptor) {
	for (int i = 0; i < _cnt; i++) {
		srand(cnt[i]);
		for (int j = 0; j < bit_of_ID;j++) {
			long_coeff[i][j] = rand() % 2;
		}
	}

	for (int i = 0; i < _cnt; i++) {
		long_raptor[i] = raptor[cnt[i]];
	}

}

//record the data to STBF
void record(int _m, int i, STBF *s, unsigned int id_num, ID myid) {

	int h1 = hash_1(_m, id_num);
	int h2 = hash_2(_m, id_num);
	int h3 = hash_3(_m, id_num);
	unsigned short ho = hash_o(id_num);
	unsigned char ho_1 = hash_o_1(id_num);
	bool myraptor = encode(myid, i);
	//no information
	if (s[i * _m + h1].flag == 0 && s[i * _m + h1].raptor == 0 && s[i * _m + h1].hash_print == 0 && s[i * _m + h1].hash_print_1 == 0) {
		s[i * _m + h1].flag = 1;
		s[i * _m + h1].raptor = myraptor;
		s[i * _m + h1].hash_print = ho;
		s[i * _m + h1].hash_print_1 = ho_1;
	}
	//collision
	else if (s[i * _m + h1].flag == 1 && (s[i * _m + h1].raptor != myraptor || s[i * _m + h1].hash_print != ho || s[i * _m + h1].hash_print_1 != ho_1)) {
		s[i * _m + h1].flag = 0;
		s[i * _m + h1].raptor = 1;
		s[i * _m + h1].hash_print = hash_print_max;
		s[i * _m + h1].hash_print_1 = hash_print_1_max;
	}

	if (s[i * _m + h2].flag == 0 && s[i * _m + h2].raptor == 0 && s[i * _m + h2].hash_print == 0 && s[i * _m + h2].hash_print_1 == 0) {
		s[i * _m + h2].flag = 1;
		s[i * _m + h2].raptor = myraptor;
		s[i * _m + h2].hash_print = ho;
		s[i * _m + h2].hash_print_1 = ho_1;
	}
	else if (s[i * _m + h2].flag == 1 && (s[i * _m + h2].raptor != myraptor || s[i * _m + h2].hash_print != ho || s[i * _m + h2].hash_print_1 != ho_1)) {
		s[i * _m + h2].flag = 0;
		s[i * _m + h2].raptor = 1;
		s[i * _m + h2].hash_print = hash_print_max;
		s[i * _m + h2].hash_print_1 = hash_print_1_max;
	}

	if (s[i * _m + h3].flag == 0 && s[i * _m + h3].raptor == 0 && s[i * _m + h3].hash_print == 0 && s[i * _m + h3].hash_print_1 == 0) {
		s[i * _m + h3].flag = 1;
		s[i * _m + h3].raptor = myraptor;
		s[i * _m + h3].hash_print = ho;
		s[i * _m + h3].hash_print_1 = ho_1;
	}
	else if (s[i * _m + h3].flag == 1 && (s[i * _m + h3].raptor != myraptor || s[i * _m + h3].hash_print != ho || s[i * _m + h3].hash_print_1 != ho_1)) {
		s[i * _m + h3].flag = 0;
		s[i * _m + h3].raptor = 1;
		s[i * _m + h3].hash_print = hash_print_max;
		s[i * _m + h3].hash_print_1 = hash_print_1_max;
	}
}


//give a group number to each cell
void form_group_col(vector<vector<int> > &G, int _m, int col_i, STBF *s) {
	for (int i = 0;i < a;i++) {
		if (s[i * _m + col_i].flag == 1 && hash_print_memory[s[i * _m + col_i].hash_print] == 0) {
			if (G.empty()) {
				vector<int> tmp;
				int tmp_p;
				tmp_p = i;
				tmp.push_back(tmp_p);
				G.push_back(tmp);
			}
			else {
				bool flag = 0;
				vector<vector<int> >::iterator iter;
				for (iter = G.begin(); iter != G.end(); iter++) {
					if (s[i * _m + col_i].hash_print == s[*iter->begin() * _m + col_i].hash_print && s[i * _m + col_i].hash_print_1 == s[*iter->begin() * _m + col_i].hash_print_1) {
						iter->push_back(i);
						flag = 1;
						break;
					}
				}
				if (flag == 0) {
					vector<int> tmp;
					int tmp_p;
					tmp_p = i;
					tmp.push_back(tmp_p);
					G.push_back(tmp);
				}
			}
		}
	}
}

//to find how many collisions 
int find_colli(int _m, int col_i, STBF *s) {
	int cnt = 0;
	for (int i = 0;i < a;i++) {
		if (s[i * _m + col_i].flag == 0 && s[i * _m + col_i].hash_print == hash_print_max && s[i * _m + col_i].raptor == 1) {
			cnt++;
		}
	}
	return cnt;
}

void record_rows(int &_cnt, int *cnt, int col_i, bool *raptor, vector<vector<int> >::iterator iter, STBF *s) {
	bool * memory_row = new bool[a];
	for (int II=0; II<a; II++) memory_row[II]=0;
	vector<int>::iterator iter1;
	for (iter1 = iter->begin(); iter1 != iter->end(); iter1++) {
		if (memory_row[*iter1] == false) {
			memory_row[*iter1] = true;
			raptor[*iter1] = s[*iter1 * m + col_i].raptor;
		}
	}
	for (int i = 0; i < a; i++) {
		if (memory_row[i] == 1) {
			cnt[_cnt] = i;
			_cnt++;
		}
	}
}


struct node {ID x; int y;} p[10000005];
map <ID,int> Value;
int cmp(node i,node j) {return i.y>j.y;}
int main() {
	int out_cnt = 0;
	int K;
	cin>>m>>K>>a;
	m=m*1024/4/a*a;
	b = 10000000 / a;
    data_of_id = new ID_input[a*b];
    Time = new int[a*b];
	s = new STBF[a*m];
	unsigned int i = 0;
	char ss[1005];
	unsigned int j = 0;
	for (i = 0; i < a * b; i++) {
	  fin.read((char *)(&data_of_id[i]), sizeof(data_of_id[i]));
	}
	fin.close();
	for (i = 0;i < a;i++) {
		id_set.clear();
		for (j = 0;j < b;j++) {
			ID input_my_id;
			memcpy(input_my_id.x, data_of_id[i * b + j].x, sizeof(input_my_id.x));
			unsigned int _data = murmur3_32(input_my_id.x, str_max_len, 1);
			if (cold_filter(_data)) {
				record(m, i, s, _data, input_my_id);
			}
			if (id_set.find(input_my_id) == id_set.end()) {
				id_set.insert(input_my_id);
				all_id_map[input_my_id] += 1;
			}
		}
	}
	map<ID, int>::iterator it;
	int CNT=0;
	int acc_cnt = 0;
	for (it = all_id_map.begin();it != all_id_map.end();it++) {
		{
		    p[++CNT].x=it->first;
		    p[CNT].y=it->second;
		}
	}
	sort(p+1,p+CNT+1,cmp);
	for (i=1; i<=K; i++) Value[p[i].x]=p[i].y;
	int col_i = 0;
	for (col_i = 0; col_i < m; col_i++) {
		int num_coli = find_colli(m, col_i, s);
		vector<vector<int> > G;	
		form_group_col(G, m, col_i, s);
		vector<vector<int> >::iterator iter;
		for (iter = G.begin();iter != G.end();iter++) {
			if (iter->size()  > g) {
				int _cnt = 0;
				int *cnt = new int[a] ;
				bool *raptor = new bool[a];
				for (int II=0; II<a; II++) cnt[II]=0,raptor[II]=0;
				record_rows(_cnt, cnt, col_i, raptor, iter, s);
				int coli[2] = { 0 };
				coli[0] = -1;
				coli[1] = -1;

				for (i = 0; i < a; i++) {
					if (s[i * m + col_i].flag == 0 && s[i * m + col_i].hash_print == hash_print_max && s[i * m + col_i].hash_print_1 == hash_print_1_max) {
						if (coli[0] == -1 || coli[1] == -1) {
							for (j = 0; j < m;j++) {
								if (s[*iter->begin() * m + col_i].hash_print == s[i * m + j].hash_print && s[*iter->begin() * m + col_i].hash_print_1 == s[i * m + j].hash_print_1) {
									cnt[_cnt] = i;
									raptor[i] = s[i * m + j].raptor;
									_cnt++;

									if (coli[0] == -1) {
										coli[0] = j;
									}
									else {
										if (coli[0] != j && coli[1] == -1) {
											coli[1] = j;
										}
									}
									break;
								}
							}
						}
						else {
							for (int coli_i = 0; coli_i < 2;coli_i++) {
								if (s[*iter->begin() * m + col_i].hash_print == s[i * m + coli[coli_i]].hash_print && s[*iter->begin() * m + col_i].hash_print_1 == s[i * m + coli[coli_i]].hash_print_1) {
									cnt[_cnt] = i;
									raptor[i] = s[i * m + coli[coli_i]].raptor;
									_cnt++;
									break;
								}
							}
						}
					}
				}

				if (_cnt > g) {
					bool *long_raptor = new bool[_cnt]();


					bool(*long_coeff)[bit_of_ID] = new bool[_cnt][bit_of_ID]();
					parameter(_cnt, cnt, raptor, long_coeff, long_raptor);

					ID myid = gauss(_cnt, long_coeff, long_raptor);
					unsigned int id_num = murmur3_32(myid.x, str_max_len, 1);
					bool test_flag = 1;
					unsigned short test_ho = hash_o(id_num);
					unsigned char test_ho_1 = hash_o_1(id_num);
					int h1 = hash_1(m, id_num);
					int h2 = hash_2(m, id_num);
					int h3 = hash_3(m, id_num);
					if (test_ho == s[*iter->begin() * m + col_i].hash_print && test_ho_1 == s[*iter->begin() * m + col_i].hash_print_1) {
						if (col_i != h1 && col_i != h2 && col_i != h3) {
							test_flag = 0;
						}
					}
					else {
						test_flag = 0;
					}

					if (test_flag) {
						out_cnt++;
						id_map.insert(make_pair(myid, _cnt));
						hash_print_memory[s[*iter->begin() * m + col_i].hash_print] = 1;
					}
					else {
					}

					delete[] long_raptor;
					delete[] long_coeff;
				}
			}
		}
	}
	CNT=0;
	for (multimap<ID, int>::iterator iter = id_map.begin(); iter != id_map.end(); iter++) {

		double x1 = all_id_map[iter->first];
		double x2 = iter->second;
		p[++CNT].x=iter->first;
		p[CNT].y=iter->second;
	}
	sort(p+1,p+CNT+1,cmp);
	cout<<CNT<<endl;
	for (i=CNT+1; i<=K; i++) p[i].x=p[0].x,p[i].y=p[0].y;
	int ans=0;
	for (i=1; i<=K; i++)
	{
	    if (Value[p[i].x]) ans++;
	    ARE+=abs(all_id_map[p[i].x]-p[i].y)/(all_id_map[p[i].x]+0.0);
	}
	double recall_rate = (double)out_cnt / (double)acc_cnt;
	printf("%d/%d=%.10f\n",ans,K,ans/(K+0.0));
	printf ("%.10f\n",ARE/K);
	return 0;
}
