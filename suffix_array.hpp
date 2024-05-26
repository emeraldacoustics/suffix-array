#include <utility>

using namespace std;

template <class T>
class SuffixArray
{
public:
	static const int sigma = 26;
	static const int maxn = 100000;
	static const int maxh = 17; // 2 ^ (maxh - 1) >= maxn

	int n;
	T s[maxn];
	int sa[maxn];
	int t0[maxn], t1[maxn];
	int c[maxn];
	int rnk[maxn];
	int lcpf[maxn];
	int rmqf[maxn][maxh];

	int idx(const T & ch)
	{
		return ch - 'a';
	}

	bool is_equal(const int * r, const int & a, const int & b, const int & l)
	{
		int la = r[a], lb = r[b], ra, rb;
		ra = a + l < n ? r[a + l] : -1;
		rb = b + l < n ? r[b + l] : -1;
		return (la == lb) && (ra == rb);
	}

	void build(const T * first, const T * last, const int & sgm = sigma)
	{
		const T * & s = first;
		n = last - first;
		memcpy(this->s, s, sizeof s[0] * n);

		int * x = t0, * y = t1, m = sgm;
		memset(c, 0, sizeof c[0] * m);
		for (int i = 0; i < n; i++)
			c[x[i] = idx(s[i])]++;
		for (int i = 1; i < m; i++)
			c[i] += c[i - 1];
		for (int i = n - 1; i >= 0; i--)
			sa[--c[x[i]]] = i;
		for (int k = 1; k <= n; k <<= 1)
		{
			int p = 0;
			for (int i = n - k; i < n; i++)
				y[p++] = i;
			for (int i = 0; i < n; i++)
			{
				if (sa[i] >= k)
					y[p++] = sa[i] - k;
			}
			memset(c, 0, sizeof c[0] * m);
			for (int i = 0; i < n; i++)
				c[x[y[i]]]++;
			for (int i = 1; i < m; i++)
				c[i] += c[i - 1];
			for (int i = n - 1; i >= 0; i--)
				sa[--c[x[y[i]]]] = y[i];
			swap(x, y);
			p = 1;
			x[sa[0]] = 0;
			for (int i = 1; i < n; i++)
				x[sa[i]] = is_equal(y, sa[i - 1], sa[i], k) ? p - 1 : p++;
			if (p >= n)
				break;
			m = p;
		}

		//get height
		for (int i = 0; i < n; i++)
			rnk[sa[i]] = i;
		lcpf[0] = 0;
		for (int i = 0, k = 0; i < n; i++)
		{
			if (rnk[i] == 0)
				continue;
			if (k > 0)
				k--;
			int u = sa[rnk[i] - 1];
			for (; i + k < n && u + k < n && s[i + k] == s[u + k]; k++);
			lcpf[rnk[i]] = k;
		}

		//get range minimum query
		for (int i = 0, k; i < n; i++)
		{
			for (k = 0; 1 << k <= i; k++)
				rmqf[i][k] = k == 0 ? lcpf[i] : min(rmqf[i][k - 1], rmqf[i - (1 << k - 1)][k - 1]);
			for (; k < maxh; k++)
				rmqf[i][k] = 0;
		}
	}

	int lcp(const int & lhs, const int & rhs)
	{
		if (lhs == rhs)
			return n - lhs;
		else
		{
			int l = rnk[lhs], r = rnk[rhs];
			if (l > r)
				swap(l, r);
			int k = 31 - __builtin_clz(r - l);
			return min(rmqf[r][k], rmqf[l + (1 << k)][k]);
		}
	}

	pair<int, int> prefix_range(const int & first, const int & last)
	{
		if (first >= last)
			return make_pair(0, n);
		else
		{
			const int x = first, l = last - first;
			int u = rnk[x];
			for (int k = maxh - 1; k >= 0; k--)
			{
				int b = 1 << k;
				if (u >= b && rmqf[u][k] >= l)
					u -= b;
			}
			int v = rnk[x];
			for (int k = maxh - 1; k >= 0; k--)
			{
				int b = 1 << k;
				if (v + b < n && rmqf[v + b][k] >= l)
					v += b;
			}
			return make_pair(u, ++v);
		}
	}
};
