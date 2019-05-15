#ifndef CALC3
#define CALC3


#include <cmath>
#include <algorithm>

namespace c3d
{

template<typename ValType, int NumRows, int NumCols>
class Mat;

template<typename ValType>
class Mat<ValType, 2, 1> {
public:

	ValType x;
	ValType y;

	const ValType& operator[](const int i) const { return (&x)[i]; }
	ValType& operator[](const int i) { return (&x)[i]; }
};

using vec2 = Mat<float, 2, 1>;
using ivec2 = Mat<int, 2, 1>;

template<typename ValType>
class Mat<ValType, 3, 1> {
public:

	ValType x;
	ValType y;
	ValType z;

	const ValType& operator[](const int i) const { return (&x)[i]; }
	ValType& operator[](const int i) { return (&x)[i]; }

};

using vec3 = Mat<float, 3, 1>;
using ivec3 = Mat<int, 3, 1>;

template<typename ValType>
class Mat<ValType, 4, 1> {
public:

	ValType x;
	ValType y;
	ValType z;
	ValType w;

	const ValType& operator[](const int i) const { return (&x)[i]; }
	ValType& operator[](const int i) { return (&x)[i]; }

};

using vec4 = Mat<float, 4, 1>;
using ivec4 = Mat<int, 4, 1>;

template<typename ValType>
class Mat<ValType, 2, 2> {
public:

	Mat<ValType, 2, 1> col_vec_0;
	Mat<ValType, 2, 1> col_vec_1;

	const auto& operator[](const int i) const { return (&col_vec_0)[i]; }
	auto& operator[](const int i) { return (&col_vec_0)[i]; }

	auto& T()
	{
		std::swap((*this)[0][1], (*this)[1][0]);
		return *this;
	}
};

using imat2 = Mat<int, 2, 2>;
using mat2 = Mat<float, 2, 2>;

template<typename ValType>
class Mat<ValType, 3, 3> {
public:

	Mat<ValType, 3, 1> col_vec_0;
	Mat<ValType, 3, 1> col_vec_1;
	Mat<ValType, 3, 1> col_vec_2;

	const auto& operator[](const int i) const { return (&col_vec_0)[i]; }
	auto& operator[](const int i) { return (&col_vec_0)[i]; }

	auto& T()
	{
		std::swap((*this)[0][1], (*this)[1][0]);
		std::swap((*this)[0][2], (*this)[2][0]);
		std::swap((*this)[1][2], (*this)[2][1]);
		return *this;
	}
};

using imat3 = Mat<int, 3, 3>;
using mat3 = Mat<float, 3, 3>;

template<typename ValType>
class Mat<ValType, 4, 4> {
public:

	Mat<ValType, 4, 1> col_vec_0;
	Mat<ValType, 4, 1> col_vec_1;
	Mat<ValType, 4, 1> col_vec_2;
	Mat<ValType, 4, 1> col_vec_3;

	const auto& operator[](const int i) const { return (&col_vec_0)[i]; }
	auto& operator[](const int i) { return (&col_vec_0)[i]; }

	auto& T()
	{
		std::swap((*this)[0][1], (*this)[1][0]);
		std::swap((*this)[0][2], (*this)[2][0]);
		std::swap((*this)[0][3], (*this)[3][0]);
		std::swap((*this)[1][2], (*this)[2][1]);
		std::swap((*this)[1][3], (*this)[3][1]);
		std::swap((*this)[2][3], (*this)[3][2]);
		return *this;
	}
};

using imat4 = Mat<int, 4, 4>;
using mat4 = Mat<float, 4, 4>;

////
// begin & end functions:
// Iterating each values in Mat in column major order.
////

template<typename ValType, int NumRows, int NumCols>
std::enable_if_t<NumCols == 1, ValType*>
	begin(Mat<ValType, NumRows, NumCols>& a)
{
	return (ValType*)(&a[0]);
}

template<typename ValType, int NumRows, int NumCols>
std::enable_if_t<NumCols != 1, ValType*>
	begin(Mat<ValType, NumRows, NumCols>& a)
{
	return (ValType*)(&a[0][0]);
}

template<typename ValType, int NumRows, int NumCols>
std::enable_if_t<NumCols == 1, const ValType*>
	begin(const Mat<ValType, NumRows, NumCols>& a)
{
	return (const ValType*)(&a[0]);
}

template<typename ValType, int NumRows, int NumCols>
std::enable_if_t<NumCols != 1, const ValType*>
	begin(const Mat<ValType, NumRows, NumCols>& a)
{
	return (const ValType*)(&a[0][0]);
}

template<typename ValType, int NumRows, int NumCols>
ValType* end(Mat<ValType, NumRows, NumCols>& a)
{
	return begin(a) + NumRows * NumCols;
}

template<typename ValType, int NumRows, int NumCols>
const ValType* end(const Mat<ValType, NumRows, NumCols>& a)
{
	return begin(a) + NumRows * NumCols;
}

////
// Element-wise operations.
////

#define ElementwiseOpsForMat(Op) \
template<typename ValType1, typename ValType2, int NumRows, int NumCols> \
auto operator##Op##=( \
Mat<ValType1, NumRows, NumCols>& lhs, \
const Mat<ValType2, NumRows, NumCols>& rhs) \
{ \
for (int i = 0; i < NumRows * NumCols; ++i) \
	* (begin(lhs) + i) Op##= *(begin(rhs) + i); \
return lhs; \
} \
\
template<typename ValType, typename ScalarType, int NumRows, int NumCols> \
std::enable_if_t< \
	std::is_same_v<float, ScalarType> || std::is_same_v<int, ScalarType>, \
	Mat<ValType, NumRows, NumCols>> \
operator##Op##=(Mat<ValType, NumRows, NumCols> & lhs, const ScalarType & rhs) \
{ \
for (int i = 0; i < NumRows * NumCols; ++i) \
	* (begin(lhs) + i) Op##= rhs; \
return lhs; \
} \
\
template<typename ValType1, typename ValType2, int NumRows, int NumCols> \
auto operator##Op( \
const Mat<ValType1, NumRows, NumCols> & lhs, \
const Mat<ValType2, NumRows, NumCols> & rhs) \
{ \
Mat<decltype(ValType1() + ValType2()), NumRows, NumCols> tmp{}; \
for (int i = 0; i < NumRows * NumCols; ++i) \
	* (begin(tmp) + i) = *(begin(lhs) + i) Op *(begin(rhs) + i); \
return tmp; \
} \
\
template<typename ValType, typename ScalarType, int NumRows, int NumCols> \
std::enable_if_t< \
std::is_same_v<float, ScalarType> || std::is_same_v<int, ScalarType>, \
Mat<decltype(ValType() + ScalarType()), NumRows, NumCols>> \
operator##Op( \
const Mat<ValType, NumRows, NumCols> & lhs, \
const ScalarType & rhs) \
{ \
Mat<decltype(ValType() + ScalarType()), NumRows, NumCols> tmp{}; \
for (int i = 0; i < NumRows * NumCols; ++i) \
	* (begin(tmp) + i) = *(begin(lhs) + i) Op rhs; \
return tmp; \
} \
\
template<typename ValType, typename ScalarType, int NumRows, int NumCols> \
std::enable_if_t< \
std::is_same_v<float, ScalarType> || std::is_same_v<int, ScalarType>, \
Mat<decltype(ValType() * ScalarType()), NumRows, NumCols>> \
operator##Op( \
const ScalarType & lhs, \
const Mat<ValType, NumRows, NumCols> & rhs) \
{ \
Mat<decltype(ValType() + ScalarType()), NumRows, NumCols> tmp{}; \
for (int i = 0; i < NumRows * NumCols; ++i) \
	* (begin(tmp) + i) = lhs Op * (begin(rhs) + i); \
return tmp; \
}

template<typename ValType1, typename ValType2, int NumRows, int NumCols>
bool operator==(
	const Mat<ValType1, NumRows, NumCols>& lhs,
	const Mat<ValType2, NumRows, NumCols>& rhs)
{
	auto lhs_iter = begin(lhs);
	auto rhs_iter = begin(rhs);
	while (lhs_iter != end(lhs))
		if (*lhs_iter != *rhs_iter)
			return false;
	return true;
}

template<typename ValType1, typename ValType2, int NumRows, int NumCols>
bool operator!=(
	const Mat<ValType1, NumRows, NumCols>& lhs,
	const Mat<ValType2, NumRows, NumCols>& rhs)
{
	return !(lhs == rhs);
}

ElementwiseOpsForMat(+);
ElementwiseOpsForMat(-);
ElementwiseOpsForMat(*);
ElementwiseOpsForMat(/);

////
// Matrix operations.
////

inline auto Det(const Mat<float, 2, 2> &m)
{
	const auto* p = begin(m);
	return p[0] * p[3] - p[1] * p[2];
}

inline auto Det(const Mat<float, 3, 3> &m)
{
	const auto* p = begin(m);
	return p[0] * (p[4] * p[8] - p[5] * p[7])
		- p[3] * (p[1] * p[8] - p[2] * p[7])
		+ p[6] * (p[1] * p[5] - p[2] * p[4]);
}

inline auto Det(const Mat<float, 4, 4> &m)
{
	const auto* p = begin(m);
	return p[10] * p[13] * p[3] * p[4] + p[0] * p[10] * p[15] * p[5]
		- p[10] * p[12] * p[3] * p[5] + p[11] * (-p[13] * p[2] * p[4]
		- p[0] * p[14] * p[5] + p[12] * p[2] * p[5] + p[0] * p[13] * p[6])
		- p[0] * p[10] * p[13] * p[7] - p[15] * p[2] * p[5] * p[8]
		+ p[14] * p[3] * p[5] * p[8] - p[13] * p[3] * p[6] * p[8]
		+ p[13] * p[2] * p[7] * p[8] + p[1] * (p[11] * p[14] * p[4]
		- p[10] * p[15] * p[4] - p[11] * p[12] * p[6] + p[10] * p[12] * p[7]
		+ p[15] * p[6] * p[8] - p[14] * p[7] * p[8])
		+ p[15] * p[2] * p[4] * p[9] - p[14] * p[3] * p[4] * p[9]
		- p[0] * p[15] * p[6] * p[9] + p[12] * p[3] * p[6] * p[9]
		+ p[0] * p[14] * p[7] * p[9] - p[12] * p[2] * p[7] * p[9];
}

inline auto Inv(const Mat<float, 2, 2>& m)
{
	const auto* p = begin(m);
	return Mat<float, 2, 2>{
		{p[3], -p[1]},
		{ -p[2], p[0] }
	} /= Det(m);
}

inline auto Inv(const Mat<float, 3, 3>& m)
{
	const auto *p = begin(m);
	return Mat<float, 3, 3>{
	{
		-p[5] * p[7] + p[4] * p[8],
			p[2] * p[7] - p[1] * p[8],
			-p[2] * p[4] + p[1] * p[5]
	},
	{
		p[5] * p[6] - p[3] * p[8],
		-p[2] * p[6] + p[0] * p[8],
		p[2] * p[3] - p[0] * p[5]
	},
	{
		-p[4] * p[6] + p[3] * p[7],
		p[1] * p[6] - p[0] * p[7],
		-p[1] * p[3] + p[0] * p[4]
	}
	} /= Det(m);
}

inline auto Inv(const Mat<float, 4, 4>& m)
{
	const auto* p = begin(m);
	return Mat<float, 4, 4>{
	{
		-p[11] * p[14] * p[5] + p[10] * p[15] * p[5] + p[11] * p[13] * p[6]
			- p[10] * p[13] * p[7] - p[15] * p[6] * p[9] + p[14] * p[7] * p[9],
			p[1] * p[11] * p[14] - p[1] * p[10] * p[15] - p[11] * p[13] * p[2]
			+ p[10] * p[13] * p[3] + p[15] * p[2] * p[9] - p[14] * p[3] * p[9],
			-p[15] * p[2] * p[5] + p[14] * p[3] * p[5] + p[1] * p[15] * p[6]
			- p[13] * p[3] * p[6] - p[1] * p[14] * p[7] + p[13] * p[2] * p[7],
			p[11] * p[2] * p[5] - p[10] * p[3] * p[5] - p[1] * p[11] * p[6]
			+ p[1] * p[10] * p[7] + p[3] * p[6] * p[9] - p[2] * p[7] * p[9]
	},
	{
		p[11] * p[14] * p[4] - p[10] * p[15] * p[4] - p[11] * p[12] * p[6]
		+ p[10] * p[12] * p[7] + p[15] * p[6] * p[8] - p[14] * p[7] * p[8],
		-p[0] * p[11] * p[14] + p[0] * p[10] * p[15] + p[11] * p[12] * p[2]
		- p[10] * p[12] * p[3] - p[15] * p[2] * p[8] + p[14] * p[3] * p[8],
		p[15] * p[2] * p[4] - p[14] * p[3] * p[4] - p[0] * p[15] * p[6]
		+ p[12] * p[3] * p[6] + p[0] * p[14] * p[7] - p[12] * p[2] * p[7],
		-p[11] * p[2] * p[4] + p[10] * p[3] * p[4] + p[0] * p[11] * p[6]
		- p[0] * p[10] * p[7] - p[3] * p[6] * p[8] + p[2] * p[7] * p[8]
	},
	{
		-p[11] * p[13] * p[4] + p[11] * p[12] * p[5] - p[15] * p[5] * p[8]
		+ p[13] * p[7] * p[8] + p[15] * p[4] * p[9] - p[12] * p[7] * p[9],
		-p[1] * p[11] * p[12] + p[0] * p[11] * p[13] + p[1] * p[15] * p[8]
		- p[13] * p[3] * p[8] - p[0] * p[15] * p[9] + p[12] * p[3] * p[9],
		-p[1] * p[15] * p[4] + p[13] * p[3] * p[4] + p[0] * p[15] * p[5]
		- p[12] * p[3] * p[5] + p[1] * p[12] * p[7] - p[0] * p[13] * p[7],
		p[1] * p[11] * p[4] - p[0] * p[11] * p[5] + p[3] * p[5] * p[8]
		- p[1] * p[7] * p[8] - p[3] * p[4] * p[9] + p[0] * p[7] * p[9]
	},
	{
		p[10] * p[13] * p[4] - p[10] * p[12] * p[5] + p[14] * p[5] * p[8]
		- p[13] * p[6] * p[8] - p[14] * p[4] * p[9] + p[12] * p[6] * p[9],
		p[1] * p[10] * p[12] - p[0] * p[10] * p[13] - p[1] * p[14] * p[8]
		+ p[13] * p[2] * p[8] + p[0] * p[14] * p[9] - p[12] * p[2] * p[9],
		p[1] * p[14] * p[4] - p[13] * p[2] * p[4] - p[0] * p[14] * p[5]
		+ p[12] * p[2] * p[5] - p[1] * p[12] * p[6] + p[0] * p[13] * p[6],
		-p[1] * p[10] * p[4] + p[0] * p[10] * p[5] - p[2] * p[5] * p[8]
		+ p[1] * p[6] * p[8] + p[2] * p[4] * p[9] - p[0] * p[6] * p[9]
	}
	} /= Det(m);
}

template<typename ValType, int NumCols, int NumRows>
ValType ValueSum(const Mat<ValType, NumCols, NumRows>& a)
{
	ValType sum{ 0 };
	for (const auto* p = begin(a); p != end(a); ++p)
		sum += *p;
	return sum;
}

// Inner-product of two vecters.
template<int NumRows, int NumCols>
std::enable_if_t<NumCols == 1, float>
	Dot(const Mat<float, NumRows, NumCols>& p,
		const Mat<float, NumRows, NumCols>& q)
{
	auto result = p * q;
	return ValueSum(result);
}

// Inner-product of matrix A, with vector v (Av).
template<int MatSize, int NumCols>
std::enable_if_t<MatSize != 1 && NumCols == 1, Mat<float, MatSize, 1>>
	Dot(const Mat<float, MatSize, MatSize>& A,
		const Mat<float, MatSize, NumCols>& v)
{
	Mat<float, MatSize, 1> result{};
	for (int i = 0; i < MatSize; ++i)
		result += A[i]*v[i];
	return result;
}

// Inner-product of two matrices, A and B.
template<int MatSize>
std::enable_if_t<MatSize != 1, Mat<float, MatSize, MatSize>>
	Dot(const Mat<float, MatSize, MatSize>& A,
		const Mat<float, MatSize, MatSize>& B)
{
	Mat<float, MatSize, MatSize> result{};
	for (int i = 0; i < MatSize; ++i)
		result[i] = Dot(A, B[i]);
	return result;
}

template<int NumRows>
auto Length(const Mat<float,NumRows, 1>& v)
{
	return std::sqrtf(Dot(v,v));
}

template<int NumRows>
auto Normalize(const Mat<float, NumRows, 1>& v)
{
	return v / Length(v);
}

inline vec3 Cross(const vec3& p, const vec3& q)
{
	return vec3{
		p.y * q.z - q.y * p.z,
		p.z * q.x - q.z * p.x,
		p.x * q.y - q.x * p.y
	};
}

////
// Quaternion w+xi+yj+zk.
////

class quat {
public:
	float w;
	float x;
	float y;
	float z;

	quat()
		: w{1}, x{ 0 }, y{ 0 }, z{ 0 }
	{}

	quat(float w, float x, float y, float z)
		: w{ w }, x{ x }, y{ y }, z{ z }
	{}

	quat(float w, const vec3& v)
		: w{ w }, x{ v.x }, y{ v.y }, z{ v.z }
	{}

	float real() const { return w; }
	void real(float w) { this->w = w; }
	vec3 imag() const { return vec3{ x, y, z }; }
	void imag(const vec3& v) { x = v.x; y = v.y; z = v.z; }
	vec3 rotate(const vec3& v)
	{
		// TODO: compare two methods.
		auto tmp1 = 2.f * Dot(imag(), v) * imag();
		auto tmp2 = (w * w - Dot(imag(), imag())) * v;
		auto tmp3 = 2.f * w * Cross(imag(), v);
		return (tmp1 + tmp2 + tmp3) / (x*x+y*y+z*z+w*w);
		// return (((*this) * quat(v, 0)) * Inv(*this)).imag();
	}

	static quat AngleAxis(float angle, vec3 axis);
};

inline quat quat::AngleAxis(float angle, vec3 axis)
{
	return quat{ cos(angle * .5f),
		Normalize(axis) * std::sinf(angle * .5f)};
}

inline float Dot(const quat& p, const quat& q)
{
	return p.x*q.x + p.y*q.y + p.z*q.z + p.w*q.w;
}

inline quat& operator+=(quat lhs, const quat& rhs)
{
	lhs.x += rhs.x;
	lhs.y += rhs.y;
	lhs.z += rhs.z;
	lhs.w += rhs.w;
	return lhs;
}

inline quat operator+(quat lhs, const quat& rhs)
{
	return lhs += rhs;
}

inline quat& operator-=(quat lhs, const quat& rhs)
{
	lhs.x -= rhs.x;
	lhs.y -= rhs.y;
	lhs.z -= rhs.z;
	lhs.w -= rhs.w;
	return lhs;
}

inline quat operator-(quat lhs, const quat& rhs)
{
	return lhs -= rhs;
}

inline quat operator*(const quat& p, float s)
{
	return quat{p.w*s,p.x*s,p.y*s,p.z*s};
}

inline quat operator*(float s, const quat& p)
{
	return p*s;
}

inline quat& operator*=(quat & lhs, const quat & rhs)
{
	const quat lhs_{ lhs };
	lhs.w = lhs_.w * rhs.w - lhs_.x * rhs.x - lhs_.y * rhs.y - lhs_.z * rhs.z;
	lhs.x = lhs_.w * rhs.x + lhs_.x * rhs.w + lhs_.y * rhs.z - lhs_.z * rhs.y;
	lhs.y = lhs_.w * rhs.y + lhs_.y * rhs.w + lhs_.z * rhs.x - lhs_.x * rhs.z;
	lhs.z = lhs_.w * rhs.z + lhs_.z * rhs.w + lhs_.x * rhs.y - lhs_.y * rhs.x;
	return lhs;
}

inline quat operator*(quat lhs, const quat & rhs)
{
	return lhs *= rhs;
}

inline quat& operator/=(quat& p, float s)
{
	p.x /= s;
	p.y /= s;
	p.z /= s;
	p.w /= s;
	return p;
}

inline quat operator/(const quat& p, float s)
{
	quat result{ p };
	return result /= s;;
}

inline quat operator/(float w, const quat& p)
{
	return quat{w/p.x,w/p.y,w/p.z,w/p.w};
}

inline quat & operator/=(quat & lhs, const quat & rhs)
{
	const quat lhs_{ lhs };

	lhs.w = lhs_.w * rhs.w + lhs_.x * rhs.x + lhs_.y * rhs.y + lhs_.z * rhs.z;
	lhs.x = lhs_.x * rhs.w - lhs_.w * rhs.x - lhs_.z * rhs.y + lhs_.y * rhs.z;
	lhs.y = lhs_.y * rhs.w + lhs_.z * rhs.x - lhs_.w * rhs.y - lhs_.x * rhs.z;
	lhs.z = lhs_.z * rhs.w - lhs_.y * rhs.x + lhs_.x * rhs.y - lhs_.w * rhs.z;

	lhs /= Dot(rhs, rhs);

	return lhs;
}

inline quat operator/(quat lhs, const quat & rhs)
{
	return lhs /= rhs;
}

inline quat Conj(const quat& q)
{
	return quat{q.w, -q.x,-q.y,-q.z};
}

inline quat Inv(const quat& q)
{
	return Conj(q) / Dot(q,q);
}

inline float Length(const quat& q)
{
	return std::sqrtf(Dot(q,q));
}

inline quat Normalize(quat q)
{
	return q / Length(q);
}

inline mat3 quat2mat3(quat q)
{
	q = Normalize(q);

	float xx = q.x*q.x, yy = q.y*q.y, zz = q.z*q.z, xz = q.x*q.z, \
		xy = q.x*q.y, yz = q.y*q.z, wx = q.w*q.x, wy = q.w*q.y, wz = q.w*q.z;

	return mat3{
		{ 1.f - 2.f * (yy + zz), 2.f * (xy - wz), 2.f * (xz + wy)},
		{ 2.f * (xy + wz), 1.f - 2.f * (xx + zz), 2.f * (yz - wx)},
		{ 2.f * (xz - wy), 2.f * (yz + wx), 1.f - 2.f * (xx + yy)}}.T();
}

inline float* begin(quat& q) { return &(q.w); }
inline const float* begin(const quat& q) { return &(q.w); }
inline float* end(quat& q) {return begin(q) + 4; }
inline const float* end(const quat& q) {return begin(q) + 4; }

}


#endif
