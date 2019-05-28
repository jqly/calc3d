
template<int NumRows, int NumCols>
auto operator%=(
    Mat<int, NumRows, NumCols>& lhs, const Mat<int, NumRows, NumCols>& rhs)
{
    for (int i = 0; i < NumRows * NumCols; ++i)
        * (begin(lhs) + i) %= *(begin(rhs) + i);
    return lhs;
}

template<int NumRows, int NumCols>
auto operator%=(Mat<int, NumRows, NumCols> & lhs, int rhs)
{
    for (int i = 0; i < NumRows * NumCols; ++i)
        * (begin(lhs) + i) %= rhs;
    return lhs;
}

template<int NumRows, int NumCols>
auto operator%(
    Mat<int, NumRows, NumCols> lhs, const Mat<int, NumRows, NumCols>& rhs)
{
    return lhs %= rhs;
}

template<int NumRows, int NumCols>
auto operator%(Mat<int, NumRows, NumCols> lhs, int rhs)
{
    return lhs %= rhs;
}

template<int NumRows, int NumCols>
auto operator%(int lhs, const Mat<int, NumRows, NumCols> & rhs)
{
    Mat<int, NumRows, NumCols> tmp{};
    for (int i = 0; i < NumRows * NumCols; ++i)
        * (begin(tmp) + i) = lhs % *(begin(rhs) + i);
    return tmp;
}
