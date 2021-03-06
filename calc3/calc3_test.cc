#define CATCH_CONFIG_MAIN
#include "catch.hpp"

#include "calc3.h"

using namespace calc;

TEST_CASE("test vec sizes") {
	static_assert(std::is_pod_v<Vec2>, "");
	static_assert(std::is_pod_v<Vec3>, "");
	static_assert(std::is_pod_v<Vec4>, "");
	static_assert(std::is_pod_v<Mat3>, "");
	static_assert( sizeof(Vec2) == 2 * 4 , "");
	static_assert( sizeof(Vec3) == 3 * 4 , "");
	static_assert( sizeof(Vec4) == 4 * 4 , "");
	static_assert( sizeof(Mat2) == 4 * 4 , "");
	static_assert( sizeof(Mat3) == 9 * 4 , "");
	static_assert( sizeof(Mat4) == 16 * 4, "" );
	static_assert( sizeof(Quat) == 4 * 4 , "");
}

TEST_CASE("test ops") {

	Mat3 m1{ {1, 2, 3},{4, 5, 6},{7, 8, 9} };
	Mat3 m2{ {5, 7, 2},{9, 0, 1},{5, 7, 2} };
	Vec3 v1{ 4, 7, 1 };
	Vec3 v2{ 2, 4, 3 };

	auto add_m1_m2 = m1 + m2;
	Mat3 np_add_m1_m2{ {6, 9, 5},
	{13, 5, 7},
	{12, 15, 11} };
	REQUIRE(ValueSum(Abs(add_m1_m2-np_add_m1_m2)) < 1e-5f);

	auto dot_m1_v1 = Dot(m1, Normalize(v1));
	Vec3 np_dot_m1_v1{ 4.80056815, 6.27766604, 7.75476393 };
	REQUIRE(ValueSum(Abs(dot_m1_v1 - np_dot_m1_v1)) < 1e-5f);

	auto div_v1 = (Dot(m1, Normalize(v1)) * v2) * (1.f / v1);
	Vec3 np_div_v1{ 2.40028,3.5872,23.2642918 };
	REQUIRE(ValueSum(Abs(div_v1 - np_div_v1)) < 1e-4f);

	auto res = v2 * (Dot(
		(m1 + m2), 
		(Dot(m1, Normalize(v1)) + Dot(m2, Normalize(v2)) * v2)) * (1.f / v1));

	Vec3 np_res{ 385.47759403,  366.73324777, 1602.88722063 };

	REQUIRE(ValueSum(Abs(res - np_res)) < .005f);

}

TEST_CASE("test Inv, T, and Cross") {
	auto m_inv = Inv(Mat4{ {-0.07841813, -0.00387151, 0.39318398, 0.89872167},
			{1.45936991, 2.66523857, 0.23746584, -2.5529272},
			{1.29317314, -0.75439287, 1.27661454, -1.21807983},
			{1.61170184, 1.37747196, 0.2069771, -0.51187045} }.T());

	Mat4 np_m_inv{ { -0.56888733,  0.62812527,  1.38584214,  0.45946347 },
					{-0.43258395,  0.35703625,  0.43362349, -0.22591418},
					{0.11042918, -0.24417286,  0.37769743, -0.15665612},
					{0.89587593, -0.09681305, -0.62826143,  0.35261243} };
	REQUIRE(ValueSum(Abs(m_inv - np_m_inv)) < .005f);

	auto m3_t_inv = Inv(Mat3{ {0.317421  ,  1.66346629,  1.06949749},
	{-1.3187909 , -0.71138327, -0.47630731},
	{-0.88199999, -1.181229  ,  0.5775008} }.T());

	Mat3 np_m3_t_inv{ {-0.36709882,  0.44563331,  0.35084597},
	   {-0.83868346,  0.42485522, -0.41189119},
	{-0.01187842, -0.47487736,  0.74213633} };

	REQUIRE(ValueSum(Abs(m3_t_inv - np_m3_t_inv)) < .005f);

	auto m2_t_inv = Inv(
		Mat2{ {-0.79521166, -1.64446396},{-0.70630628, -0.54255584} }.T());
	Mat2 np_m2_t_inv{ {0.7431778 , -0.96747857}, {-2.2525407 ,  1.08925867} };
	REQUIRE(ValueSum(Abs(m2_t_inv - np_m2_t_inv)) < .005f);

    Vec3 p{1,2,3}, q{4,5,6};
    Vec3 cross_pq = Cross(p,q);
    Vec3 np_cross_pq{-3,6,-3};
    REQUIRE(ValueSum(Abs(cross_pq-np_cross_pq)) < 1e-4);
}

TEST_CASE("test Quat") {

	Quat q{ 1, 2, 3, 4 };
	REQUIRE(std::abs(Length(q)-std::sqrtf(30)) < 1e-5);

	Mat3 qrot{
		{-2./3., 2./3., 1./3.},
	    {2./15., -1./3., 14./15.},
	    {11./15., 2./3., 2./15.}
	};
	
	REQUIRE(ValueSum(Abs(qrot - Quat2Mat3(Normalize(q)))) < .005f);

	Quat nq = Quat::AngleAxis(
		2.f*std::acosf(1.f/std::sqrtf(30.f)), 
		Normalize(Vec3{
			std::sqrtf(2.f/15.f),
			std::sqrtf(3.f/10.f),
			2.f*std::sqrtf(2.f/15.f)}));
	REQUIRE( Length(Normalize(q) - nq) < 1e-4);

	Quat a{ 2,0,-6,3 }, b{ 1,3,-2,2 };
	Quat mul_a_b = a * b;
	REQUIRE(Length(a * b - mul_a_b) < 1e-4);

	Quat q1 = Quat::AngleAxis(3.14159f/4.f, Vec3{0,0,1});
	REQUIRE( std::abs(Length(q1) - 1.f) < 1e-6);

	Vec3 v{1,1,0};

	Vec3 rv1 = QuatRotate(q1,v);
	Vec3 rv2 = Dot(Quat2Mat3(q1),v);
	Vec3 rv3 = (q1*Quat(0,v)*Inv(q1)).Im();

	REQUIRE( Length(rv1-rv2) < 1e-4 );
	REQUIRE( Length(rv3-rv2) < 1e-4 );

}

TEST_CASE("test transformations") {
	Mat3 ma_rotated{
		{0.869853, 0.434061, -0.234407}, 
		{-0.411094, 0.900476, 0.141933}, 
		{0.272686, -0.0270977, 0.961722}};
	Mat3 rotated = RotationTransform(
		Deg2Rad(30), Normalize(Vec3({1.f,3.f,5.f})));
	
	Mat3 quat_rotated = Quat2Mat3(
		Quat::AngleAxis(Deg2Rad(30), Normalize(Vec3({1.f,3.f,5.f}))));

	REQUIRE(ValueSum(Abs(rotated - ma_rotated)) < .005f);
	REQUIRE(ValueSum(Abs(rotated - quat_rotated)) < .005f);

	Mat4 lookat = LookAt(
		{0.044106,-0.003575,6.567765}, 
		{0.044106,-0.003575,0.066791}, 
		{0.000000,1.000000,0.000000});
	Mat4 xy_lookat = Mat4{{1.000000,-0.000000,0.000000,-0.044106},
		 {0.000000,1.000000,0.000000,0.003575},
		 {-0.000000,-0.000000,1.000000,-6.567765},
		 {0.000000,0.000000,0.000000,1.000000}}.T();
	REQUIRE(ValueSum(Abs(xy_lookat - lookat)) < .005f);

	Mat4 proj = ProjectiveTransform(1.0472, 1, 3.25049, 9.75146);
	Mat4 xy_proj = Mat4{{1.732051,0.000000,0.000000,0.000000},
	{0.000000,1.732051,0.000000,0.000000},
	{0.000000,0.000000,-2.000000,-9.751461},
	{0.000000,0.000000,-1.000000,0.000000}}.T();
	REQUIRE(ValueSum(Abs(proj - xy_proj)) < .001f); 

	Mat4 ortho = OrthographicTransform(
		-2.62256,2.62256,-3.08495,3.08495,0,6.50033);
	Mat4 xy_ortho = Mat4{{0.381306,0.000000,0.000000,-0.000000},
	{0.000000,0.324154,0.000000,-0.000000},
	{0.000000,0.000000,-0.307676,-1.000000},
	{0.000000,0.000000,0.000000,1.000000}}.T();
	REQUIRE(ValueSum(Abs(ortho - xy_ortho)) < .005f); 

}

TEST_CASE("test iVec mod") {
	iMat3 m1 = iMat3{
		{21, 1, 17},
		{53, 67, 80},
		{68, 64, 13}}.T();
	iMat3 m2= iMat3{
		{16,  0, 90},
		{18, 38, 28},
		{91, 42, 90}}.T();
	iMat3 np_result= iMat3{
		{16,  0,  5},
		{18, 38, 28},
		{23, 42, 12}}.T();
	int v1 = 30, v2 = 3;
	auto result = m2%m1;
	REQUIRE( 
		std::memcmp(begin(result),begin(np_result), sizeof(float)*9) == 0 );

	auto result2 = m2%20;
	auto np_result2 = iMat3{{16,  0, 10},
       {18, 18,  8},
	{11,  2, 10}}.T();
	REQUIRE( 
		std::memcmp(begin(result2),begin(np_result2), sizeof(float)*9) == 0 );

		auto result3 = 65%m1;
	auto np_result3 = iMat3{{2,  0, 14},
       {12, 65, 65},
	{65,  1,  0}}.T();
	REQUIRE( 
		std::memcmp(begin(result3),begin(np_result3), sizeof(float)*9) == 0 );
}

TEST_CASE("test component-wise ops & string-ify") {
	iMat4 m1{
		{  8,  32,  98, -65},
		{ 33,  86,  59, -59},
		{ 54,  29,  57, -83},
		{-97, -78,  -4,  21}};
	iMat4 m2{
		{ -8,  54,  -9, -17},
		{ 26,   8,  50,  77},
		{ 74, -26,  43,  41},
		{-70,  85,  19,  18}};
	iMat4 result_minimum = ElementWiseMin(m1,m2);
	iMat4 np_result_minimum{
		{ -8,  32,  -9, -65},
		{ 26,   8,  50, -59},
		{ 54, -26,  43, -83},
		{-97, -78,  -4,  18}};

	REQUIRE(ValueSum(Abs(result_minimum - np_result_minimum)) < .005f);
	REQUIRE( std::memcmp(MatrixForm(result_minimum).data(), 
		MatrixForm(np_result_minimum).data(), sizeof(iMat4)) == 0 );
	
	Quat q{1,2,3,4};
	std::string q_result = "1.000000+2.000000i+3.000000j+4.000000k";
	REQUIRE( std::memcmp(QuaternionForm(q).data(), 
		q_result.data(), sizeof(q_result)) == 0 );

}

TEST_CASE("test cast") {
	Mat4 m4{{  8,  32,  98, -65},
		{ 33,  86,  59, -59},
		{ 54,  29,  57, -83},
		{-97, -78,  -4,  21}};
	auto m3 = Cast<Mat3>(m4);
	Mat3 m3_{{  8,  32,  98},
		{ 33,  86,  59},
		{ 54,  29,  57}};
	REQUIRE( std::memcmp(&m3, &m3_, sizeof(m3)) == 0 );

	auto v4 = Cast<Vec4>(m4);
	Vec4 v4_{8,  32,  98, -65};
	REQUIRE( std::memcmp(&v4, &v4_, sizeof(v4)) == 0 );

	auto v2 = Cast<Vec2>(v4);
	Vec4 v2_{8,  32};
	REQUIRE( std::memcmp(&v2, &v2_, sizeof(v2)) == 0 );

	auto m4_diag = Diagonal<Mat4>(1);
	Mat4 m4_diag_{{1,0,0,0},{0,1,0,0},{0,0,1,0},{0,0,0,1}};
	REQUIRE( std::memcmp( \
		&m4_diag, &m4_diag_, sizeof(m4_diag)) == 0 );
}
