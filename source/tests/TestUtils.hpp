#ifndef TEST_UTILS_HPP_
#define TEST_UTILS_HPP_

inline bool equal(const Vector3d& left, const Vector3d& right)
{
    return Approx(left.x()).margin(1e-12) == right.x()
           && Approx(left.y()).margin(1e-12) == right.y()
           && Approx(left.z()).margin(1e-12) == right.z();
}

#endif
