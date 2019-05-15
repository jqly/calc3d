#define CATCH_CONFIG_MAIN
#include "catch.hpp"

#include "calc3.h"

using namespace c3d;

TEST_CASE("test vec sizes") {
	static_assert(std::is_pod_v<vec2>, "");
	static_assert(std::is_pod_v<vec3>, "");
	static_assert(std::is_pod_v<vec4>, "");
	static_assert(std::is_pod_v<mat3>, "");
	static_assert(sizeof(vec2) == 2 * 4, "");
	static_assert(sizeof(vec3) == 3 * 4, "");
	static_assert(sizeof(vec4) == 4 * 4, "");
	static_assert(sizeof(mat2) == 4 * 4, "");
	static_assert(sizeof(mat3) == 9 * 4, "");
	static_assert(sizeof(mat4) == 16 * 4, "");
}

TEST_CASE("test ops") {

	mat3 m1{ {1, 2, 3},{4, 5, 6},{7, 8, 9} };
	mat3 m2{ {5, 7, 2},{9, 0, 1},{5, 7, 2} };
	vec3 v1{ 4, 7, 1 };
	vec3 v2{ 2, 4, 3 };

	auto add_m1_m2 = m1 + m2;
	mat3 np_add_m1_m2{ {6, 9, 5},
	{13, 5, 7},
	{12, 15, 11} };
	REQUIRE(ValueSum(add_m1_m2-np_add_m1_m2) < 1e-5f);

	auto dot_m1_v1 = Dot(m1, Normalize(v1));
	vec3 np_dot_m1_v1{ 4.80056815, 6.27766604, 7.75476393 };
	REQUIRE(ValueSum(dot_m1_v1 - np_dot_m1_v1) < 1e-5f);

	auto div_v1 = (Dot(m1, Normalize(v1)) * v2) * (1.f / v1);
	vec3 np_div_v1{ 2.40028,3.5872,23.2642918 };
	REQUIRE(ValueSum(div_v1 - np_div_v1) < 1e-4f);

	auto res = v2 * (Dot((m1 + m2), (Dot(m1, Normalize(v1)) + Dot(m2, Normalize(v2)) * v2)) * (1.f / v1));

	vec3 np_res{ 385.47759403,  366.73324777, 1602.88722063 };

	REQUIRE(ValueSum(res - np_res) < .005f);

}

TEST_CASE("test Inv, T, and Cross") {
	auto m_inv = Inv(mat4{ {-0.07841813, -0.00387151, 0.39318398, 0.89872167},
			{1.45936991, 2.66523857, 0.23746584, -2.5529272},
			{1.29317314, -0.75439287, 1.27661454, -1.21807983},
			{1.61170184, 1.37747196, 0.2069771, -0.51187045} }.T());

	mat4 np_m_inv{ { -0.56888733,  0.62812527,  1.38584214,  0.45946347 },
					{-0.43258395,  0.35703625,  0.43362349, -0.22591418},
					{0.11042918, -0.24417286,  0.37769743, -0.15665612},
					{0.89587593, -0.09681305, -0.62826143,  0.35261243} };
	REQUIRE(ValueSum(m_inv - np_m_inv) < .005f);

	auto m3_t_inv = Inv(mat3{ {0.317421  ,  1.66346629,  1.06949749},
	{-1.3187909 , -0.71138327, -0.47630731},
	{-0.88199999, -1.181229  ,  0.5775008} }.T());

	mat3 np_m3_t_inv{ {-0.36709882,  0.44563331,  0.35084597},
	   {-0.83868346,  0.42485522, -0.41189119},
	{-0.01187842, -0.47487736,  0.74213633} };

	REQUIRE(ValueSum(m3_t_inv - np_m3_t_inv) < .005f);

	auto m2_t_inv = Inv(mat2{ {-0.79521166, -1.64446396},{-0.70630628, -0.54255584} }.T());
	mat2 np_m2_t_inv{ {0.7431778 , -0.96747857}, {-2.2525407 ,  1.08925867} };
	REQUIRE(ValueSum(m2_t_inv - np_m2_t_inv) < .005f);

    vec3 p{1,2,3}, q{4,5,6};
    vec3 cross_pq = Cross(p,q);
    vec3 np_cross_pq{-3,6,-3};
    REQUIRE(ValueSum(cross_pq-np_cross_pq) < 1e-4);
}

TEST_CASE("test quat") {

	quat q{ 1, 2, 3, 4 };
	REQUIRE(std::abs(Length(q)-std::sqrtf(30)) < 1e-5);

	mat3 qrot{
		{-2./3., 2./3., 1./3.},
	    {2./15., -1./3., 14./15.},
	    {11./15., 2./3., 2./15.}
	};
	
	REQUIRE(ValueSum(qrot - quat2mat3(q)) < .005f);

	quat nq = quat::AngleAxis(
		2.f*std::acosf(1.f/std::sqrtf(30.f)), 
		vec3{
			std::sqrtf(2.f/15.f),
			std::sqrtf(3.f/10.f),
			2.f*std::sqrtf(2.f/15.f)});
	REQUIRE( Length(Normalize(q) - nq) < 1e-4);

	quat a{ 2,0,-6,3 }, b{ 1,3,-2,2 };
	quat mul_a_b = a * b;
	REQUIRE(Length(a * b - mul_a_b) < 1e-4);

}