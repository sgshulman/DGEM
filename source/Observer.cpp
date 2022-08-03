#include <cstdio>
#include "Observer.hpp"

namespace
{
    void formatAngle(char *str, int const length, double const angleDeg)
    {
        long const angleMilliDeg = std::lround(degrees(angleDeg) * 1000);
        long const degrees = angleMilliDeg / 1000 % 360;
        long const milliDegrees = angleMilliDeg % 1000;

        if (milliDegrees == 0)
        {
            snprintf(str, length, "%2.2li", degrees);
        } else
        {
            snprintf(str, length, "%2.2li-%3.3li", degrees, milliDegrees);
        }
    }
}

// pictures
Pictures::Pictures(std::uint32_t nx, std::uint32_t ny)
    : nx_{ nx }
    , ny_{ ny }
{
    f_ = new double[ nx*ny ]();
    q_ = new double[ nx*ny ]();
    u_ = new double[ nx*ny ]();
}

Pictures::~Pictures()
{
    delete[] f_;
    delete[] q_;
    delete[] u_;
}

Pictures::Pictures(Pictures&& other) noexcept
    : nx_{ other.nx_ }
    , ny_{ other.ny_ }
    , f_{ other.f_ }
    , q_{ other.q_ }
    , u_{ other.u_ }
{
    other.f_ = nullptr;
    other.q_ = nullptr;
    other.u_ = nullptr;
}

Pictures& Pictures::operator=(Pictures&& other) noexcept
{
    std::swap(f_, other.f_);
    std::swap(q_, other.q_);
    std::swap(u_, other.u_);
    nx_ = other.nx_;
    ny_ = other.ny_;
    return *this;
}

// place photon on the images
void Pictures::bin(Photon const& ph, int64_t const xl, int64_t const yl, int64_t const id, double const weight)
{
    // place weighted photon into image location
    if((id == -1 || static_cast<std::uint64_t>(id) == ph.nscat()) && (xl>=0) && (yl >= 0) && ((std::uint32_t)xl < nx_) && ((std::uint32_t)yl < ny_) )
    {
        f_[xl+yl*nx_] += weight * ph.weight() * ph.fi();
        q_[xl+yl*nx_] += weight * ph.weight() * ph.fq();
        u_[xl+yl*nx_] += weight * ph.weight()*ph.fu();
    }
}

void Pictures::normalize(std::uint64_t const numPhotons)
{
    for (std::uint64_t i=0; i!=nx_*ny_; ++i)
    {
        f_[i] /= numPhotons;
        q_[i] /= numPhotons;
        u_[i] /= numPhotons;
    }
}

void Pictures::write(double const phi, double const theta, int const key) const
{
    int const FILENAME_LENGTH{ 40 };
    int const ANGLE_LENGTH{ 10 };

    char fname[FILENAME_LENGTH];
    char phiStr[ANGLE_LENGTH];
    char thetaStr[ANGLE_LENGTH];

    formatAngle(phiStr, ANGLE_LENGTH, phi);
    formatAngle(thetaStr, ANGLE_LENGTH, theta);

    snprintf(fname, FILENAME_LENGTH, "fimage%s_%s_%2.2i.dat", phiStr, thetaStr, key);
    std::ofstream f(fname);
    f.precision(14);

    for (std::uint64_t y=0; y != ny_; ++y)
    {
        for (std::uint64_t x=0; x!=nx_; ++x)
        {
            f << f_[x + y * nx_] << "\t";
        }
        f << "\n";
    }
    f.close();

    snprintf(fname, FILENAME_LENGTH, "qimage%s_%s_%2.2i.dat", phiStr, thetaStr, key);
    std::ofstream q(fname);
    q.precision(14);

    for (std::uint64_t y=0; y != ny_; ++y)
    {
        for (std::uint64_t x=0; x!=nx_; ++x)
        {
            q << q_[x + y * nx_] << "\t";
        }
        q << "\n";
    }
    q.close();

    snprintf(fname, FILENAME_LENGTH, "uimage%s_%s_%2.2i.dat", phiStr, thetaStr, key);
    std::ofstream u(fname);
    u.precision(14);

    for (std::uint64_t y=0; y != ny_; ++y)
    {
        for (std::uint64_t x=0; x!=nx_; ++x)
        {
            u << u_[x + y * nx_] << "\t";
        }
        u << "\n";
    }
    u.close();
}

void Pictures::sum(std::ofstream& file)
{
    double fsum=0, qsum=0, usum=0;
    for (std::uint64_t x=0; x!=nx_; ++x)
    {
        for (std::uint64_t y=0; y != ny_; ++y)
        {
            fsum += f_[y + x * ny_];
            qsum += q_[y + x * ny_];
            usum += u_[y + x * ny_];
        }
    }
    file << "\tF= " << fsum << "\tQ= " << qsum << "\tU= " << usum
         << "\tp= " << std::sqrt(qsum*qsum + usum*usum)/fsum
         << "\tphi= " << 90 * std::atan2(usum, qsum)/PI;
}


Observer::Observer(double const phi, double const theta, double const rimage, std::uint32_t const Nx, std::uint32_t const Ny)
    : result_(Nx, Ny)
    , result0_(Nx, Ny)
    , result1_(Nx, Ny)
    , result2_(Nx, Ny)
    , direction_{phi, theta}
    , nx_(Nx)
    , ny_(Ny)
{
    rimage_ = rimage;
    theta_ = theta;

    cosp_ = std::cos(phi);
    sinp_ = std::sin(phi);
}

void Observer::normalize(std::uint64_t const numPhotons)
{
    result_.normalize(numPhotons);
    result0_.normalize(numPhotons);
    result1_.normalize(numPhotons);
    result2_.normalize(numPhotons);
}

void Observer::writeToMapFiles(bool const fWriteSingleAndDoubleScatterings, std::uint32_t const numberOfScatterings)
{
    result_.write(direction_.phi(), theta_, 0);

    if (fWriteSingleAndDoubleScatterings)
    {
        result1_.write(direction_.phi(), theta_, 1);

        if (numberOfScatterings >= 2)
        {
            result2_.write(direction_.phi(), theta_, 2);
        }
    }
}

void Observer::write(std::ofstream& file)
{
    file << "phi= " << degrees(direction_.phi()) << "\ttheta= " << degrees(theta_);
    result_.sum(file);
    result0_.sum(file);
    result1_.sum(file);
    result2_.sum(file);
    file << "\n";
}

bool Observer::inFov(const Photon &photon) const
{
    double const yimage = rimage_ + photon.pos().z() * direction_.sinTheta() -
                          direction_.cosTheta() * (photon.pos().y()*sinp_ + photon.pos().x()*cosp_);

    double const ximage = rimage_ + photon.pos().y()*cosp_ - photon.pos().x()*sinp_;

    auto const xl = static_cast<int64_t>(nx_ * ximage / (2.0 * rimage_));
    auto const yl = static_cast<int64_t>(ny_ * yimage / (2.0 * rimage_));

    return (xl>=0) && (yl >= 0) && ((std::uint32_t)xl < nx_) && ((std::uint32_t)yl < ny_);
}

void Observer::bin(Photon const& photon)
{
    double const yimage = imageY(photon.pos());
    double const ximage = imageX(photon.pos());

    auto const xl = static_cast<int64_t>(nx_ * ximage / (2.0 * rimage_));
    auto const yl = static_cast<int64_t>(ny_ * yimage / (2.0 * rimage_));

    double const eps = std::numeric_limits<float>::epsilon();

    bin(
        photon,
        xl,
        yl,
        photon.nscat() != 0 && std::abs(ximage - xl * (2.0 * rimage_) / nx_) < eps,
        photon.nscat() != 0 && std::abs(yimage - yl * (2.0 * rimage_) / ny_) < eps,
        1.0);
}

// TODO: Add tests for this method
void Observer::bin(Photon const& photon, const Vector3d &pos1, const Vector3d &pos2)
{
    double const yimage1 = imageY(pos1);
    double const ximage1 = imageX(pos1);
    auto const xl1 = static_cast<int64_t>(nx_ * ximage1 / (2.0 * rimage_));
    auto const yl1 = static_cast<int64_t>(ny_ * yimage1 / (2.0 * rimage_));

    double const yimage2 = imageY(pos2);
    double const ximage2 = imageX(pos2);
    auto const xl2 = static_cast<int64_t>(nx_ * ximage2 / (2.0 * rimage_));
    auto const yl2 = static_cast<int64_t>(ny_ * yimage2 / (2.0 * rimage_));

    double const eps = std::numeric_limits<float>::epsilon();

    if (xl1 == xl2 && yl1 == yl2)
    {
        bin(
            photon,
            xl2,
            yl2,
            std::abs(ximage2 - ximage1) < eps && std::abs(ximage2 - xl2 * (2.0 * rimage_) / nx_) < eps,
            std::abs(yimage2 - yimage1) < eps && std::abs(yimage2 - yl2 * (2.0 * rimage_) / ny_) < eps,
            1.0);

        return;
    }

    double xImageMin = std::min(ximage1, ximage2);
    double xImageMax = std::max(ximage1, ximage2);
    int64_t borderX = xl1;
    int64_t lastBorderX = xl2;

    if (ximage1 > ximage2)
    {
        xImageMin = ximage2;
        xImageMax = ximage1;
        borderX = xl2;
        lastBorderX = xl1;
    }

    double yImageMin = yimage1;
    double yImageMax = yimage2;
    int64_t borderY = yl1;
    int64_t lastBorderY = yl2;

    if (yimage1 > yimage2)
    {
        yImageMin = yimage2;
        yImageMax = yimage1;
        borderY = yl2;
        lastBorderY = yl1;
    }

    double const dx = xImageMax - xImageMin;
    double const dy = yImageMax - yImageMin;
    double const totalW = std::sqrt(dx*dx + dy*dy);

    double x = xImageMin;
    double y = yImageMin;

    while (borderX <= lastBorderX && borderY <= lastBorderY)
    {
        double const xborder = std::min(static_cast<double>(borderX+1) * (2.0 * rimage_) / nx_, xImageMax);
        double const yborder = std::min(static_cast<double>(borderY+1) * (2.0 * rimage_) / ny_, yImageMax);

        double const xt = dx > 0 ? (xborder - x) / dx : std::numeric_limits<double>::infinity();
        double const yt = dy > 0 ? (yborder - y) / dy : std::numeric_limits<double>::infinity();

        if (xt < yt)
        {
            double yNew = y + xt * dy;
            double w = std::sqrt((yNew - y)*(yNew-y) + (xborder - x)*(xborder - x));

            bin(
                photon,
                borderX,
                borderY,
                false,
                dy < eps && std::fabs(y - borderY * (2.0 * rimage_) / ny_) < eps,
                w / totalW);

            ++borderX;
            x = xborder;
            y = yNew;
        } else {
            double xNew = x + yt * dx;
            double w = std::sqrt((xNew - x)*(xNew-x) + (yborder - y)*(yborder - y));

            bin(
                photon,
                borderX,
                borderY,
                dx < eps && std::fabs(x - borderX * (2.0 * rimage_) / nx_) < eps,
                false,
                w / totalW);

            ++borderY;
            x = xNew;
            y = yborder;
        }
    }
}


void Observer::binHex(Photon const& photon, Vector3d const& pos1, Vector3d const& pos2)
{
    constexpr double a = 1.0 / std::sqrt(7);

    constexpr Vector3d shifts[7] = {{0., 0., 0.}, {1.5*a, std::sqrt(3.)/2.*a, 0}, {1.5*a, -std::sqrt(3.)/2.*a, 0},
        {0, std::sqrt(3.)*a, 0}, {0, -std::sqrt(3.)*a, 0},
        {-1.5*a, std::sqrt(3.)/2.*a, 0}, {-1.5*a, -std::sqrt(3.)/2.*a, 0}};

    Vector3d const v = pos2 - pos1;
    double const r = v.norm() / 2.;
    Vector3d const x = Vector3d(-v.z(), v.x() * v.y() > 0 ? -v.z() : v.z(), v.x() * v.y() > 0 ? v.x() + v.y() : v.x() - v.y()).normalized();
    Vector3d const y = vectorProduct(v, x).normalized();

    for (int i=0; i!=7; ++i)
    {
        Vector3d const shift = r * (shifts[i].x() * x + shifts[i].y() * y);
        double const yimage1 = imageY(pos1 + shift);
        double const ximage1 = imageX(pos1 + shift);
        auto const xl1 = static_cast<int64_t>(nx_ * ximage1 / (2.0 * rimage_));
        auto const yl1 = static_cast<int64_t>(ny_ * yimage1 / (2.0 * rimage_));

        double const yimage2 = imageY(pos2 + shift);
        double const ximage2 = imageX(pos2 + shift);
        auto const xl2 = static_cast<int64_t>(nx_ * ximage2 / (2.0 * rimage_));
        auto const yl2 = static_cast<int64_t>(ny_ * yimage2 / (2.0 * rimage_));

        double const eps = std::numeric_limits<float>::epsilon();

        if (xl1 == xl2 && yl1 == yl2)
        {
            bin(
                photon,
                xl2,
                yl2,
                std::abs(ximage2 - ximage1) < eps && std::abs(ximage2 - xl2 * (2.0 * rimage_) / nx_) < eps,
                std::abs(yimage2 - yimage1) < eps && std::abs(yimage2 - yl2 * (2.0 * rimage_) / ny_) < eps,
                1.0/7.0);
        } else {
            double xImageMin = std::min(ximage1, ximage2);
            double xImageMax = std::max(ximage1, ximage2);
            int64_t borderX = xl1;
            int64_t lastBorderX = xl2;

            if (ximage1 > ximage2)
            {
                xImageMin = ximage2;
                xImageMax = ximage1;
                borderX = xl2;
                lastBorderX = xl1;
            }

            double yImageMin = yimage1;
            double yImageMax = yimage2;
            int64_t borderY = yl1;
            int64_t lastBorderY = yl2;

            if (yimage1 > yimage2)
            {
                yImageMin = yimage2;
                yImageMax = yimage1;
                borderY = yl2;
                lastBorderY = yl1;
            }

            double const dx = xImageMax - xImageMin;
            double const dy = yImageMax - yImageMin;
            double const totalW = std::sqrt(dx*dx + dy*dy);

            double x = xImageMin;
            double y = yImageMin;

            while (borderX <= lastBorderX && borderY <= lastBorderY)
            {
                double const xborder = std::min(static_cast<double>(borderX+1) * (2.0 * rimage_) / nx_, xImageMax);
                double const yborder = std::min(static_cast<double>(borderY+1) * (2.0 * rimage_) / ny_, yImageMax);

                double const xt = dx > 0 ? (xborder - x) / dx : std::numeric_limits<double>::infinity();
                double const yt = dy > 0 ? (yborder - y) / dy : std::numeric_limits<double>::infinity();

                if (xt < yt)
                {
                    double yNew = y + xt * dy;
                    double w = std::sqrt((yNew - y)*(yNew-y) + (xborder - x)*(xborder - x));

                    bin(
                        photon,
                        borderX,
                        borderY,
                        false,
                        dy < eps && std::fabs(y - borderY * (2.0 * rimage_) / ny_) < eps,
                        w / totalW / 7.0);

                    ++borderX;
                    x = xborder;
                    y = yNew;
                } else {
                    double xNew = x + yt * dx;
                    double w = std::sqrt((xNew - x)*(xNew-x) + (yborder - y)*(yborder - y));

                    bin(
                        photon,
                        borderX,
                        borderY,
                        dx < eps && std::fabs(x - borderX * (2.0 * rimage_) / nx_) < eps,
                        false,
                        w / totalW / 7.0);

                    ++borderY;
                    x = xNew;
                    y = yborder;
                }
            }
        }
    }
}


inline double Observer::imageX(Vector3d const &position) const
{
    return rimage_ + position.y() * cosp_ - position.x() * sinp_;
}

inline double Observer::imageY(Vector3d const &position) const
{
    return rimage_ + position.z() * direction_.sinTheta() - direction_.cosTheta() * (position.y()*sinp_ + position.x()*cosp_);
}

inline void Observer::bin(Photon const& photon, int64_t const x, int64_t const y, double const weight)
{
    result_.bin(photon, x, y, -1, weight);
    result0_.bin(photon, x, y, 0, weight);
    result1_.bin(photon, x, y, 1, weight);
    result2_.bin(photon, x, y, 2, weight);
}

inline void Observer::bin(const Photon &photon, int64_t const x, int64_t const y, bool const onXBorder, bool const onYBorder, double const weight)
{
    if (onXBorder && x > 0)
    {
        bin(photon, x, y, 0.5 * weight);
        bin(photon, x-1, y, 0.5 * weight);
    } else if (onYBorder && y > 0) {
        bin(photon, x, y, 0.5 * weight);
        bin(photon, x, y-1, 0.5 * weight);
    } else {
        bin(photon, x, y, weight);
    }
}
