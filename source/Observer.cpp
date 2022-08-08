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
    f_ = new double[ 3*nx*ny ]();
}

Pictures::~Pictures()
{
    delete[] f_;
}

Pictures::Pictures(Pictures&& other) noexcept
    : nx_{ other.nx_ }
    , ny_{ other.ny_ }
    , f_{ other.f_ }
{
    other.f_ = nullptr;
}

Pictures& Pictures::operator=(Pictures&& other) noexcept
{
    std::swap(f_, other.f_);
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
        f_[3*(xl+yl*nx_)] += weight * ph.weight() * ph.fi();
        f_[3*(xl+yl*nx_)+1] += weight * ph.weight() * ph.fq();
        f_[3*(xl+yl*nx_)+2] += weight * ph.weight() * ph.fu();
    }
}

void Pictures::normalize(std::uint64_t const numPhotons)
{
    for (std::uint64_t i=0; i!=3*nx_*ny_; ++i)
    {
        f_[i] /= numPhotons;
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
            f << f_[3*(x + y * nx_)] << "\t";
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
            q << f_[3*(x + y * nx_)+1] << "\t";
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
            u << f_[3*(x + y * nx_)+2] << "\t";
        }
        u << "\n";
    }
    u.close();
}

void Pictures::sum(std::ostream& stream)
{
    double fsum=0, qsum=0, usum=0;
    for (std::uint64_t x=0; x!=nx_; ++x)
    {
        for (std::uint64_t y=0; y != ny_; ++y)
        {
            fsum += f_[3*(y + x * ny_)];
            qsum += f_[3*(y + x * ny_)+1];
            usum += f_[3*(y + x * ny_)+2];
        }
    }
    stream << "\tF= " << fsum << "\tQ= " << qsum << "\tU= " << usum
         << "\tp= " << std::sqrt(qsum*qsum + usum*usum)/fsum
         << "\tphi= " << 90 * std::atan2(usum, qsum)/PI;
}


Observer::Observer(double const phi, double const theta, double const rImage, double const rMask, std::uint32_t const Nx, std::uint32_t const Ny)
    : result_(Nx, Ny)
    , result0_(Nx, Ny)
    , result1_(Nx, Ny)
    , result2_(Nx, Ny)
    , direction_{phi, theta}
    , nx_{ Nx }
    , ny_{ Ny }
    , rImage_{ rImage }
    , rImageRev_{ 1.0 / (2.0 * rImage) }
    , pixelSize_{ (2.0 * rImage_) / nx_, (2.0 * rImage_) / ny_ }
    , rMask2_{ rMask * rMask }
    , theta_{ theta }
    , cosp_{ std::cos(phi) }
    , sinp_{ std::sin(phi) }
{}


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

void Observer::write(std::ostream& stream)
{
    stream << "phi= " << degrees(direction_.phi()) << "\ttheta= " << degrees(theta_);
    result_.sum(stream);
    result0_.sum(stream);
    result1_.sum(stream);
    result2_.sum(stream);
    stream << "\n";
}

bool Observer::inFov(const Vector3d &pos) const
{
    Vector2d const imagePos = project(pos);

    double const r2 = imagePos.norm2();

    auto const xl = static_cast<int64_t>(nx_ * (rImage_ + imagePos.x()) * rImageRev_);
    auto const yl = static_cast<int64_t>(ny_ * (rImage_ + imagePos.y()) * rImageRev_);

    return (xl>=0) && (yl >= 0) && ((std::uint32_t)xl < nx_) && ((std::uint32_t)yl < ny_) && r2 >= rMask2_;
}


void Observer::bin(Photon const& photon)
{
    Vector2d const imagePos = project(photon.pos());
    double const r2 = imagePos.norm2();

    if (r2 >= rMask2_)
    {
        Vector2d const imagePosShifted = imagePos + Vector2d{rImage_, rImage_};
        double const x = nx_ * imagePosShifted.x() * rImageRev_;
        double const y = ny_ * imagePosShifted.y() * rImageRev_;
        auto const xl = static_cast<std::int64_t>(x) - (x < 0.0);
        auto const yl = static_cast<std::int64_t>(y) - (y < 0.0);

        double const eps = std::numeric_limits<float>::epsilon();

        binPoint(
            photon,
            xl,
            yl,
            photon.nscat() != 0 && std::abs(imagePosShifted.x() - xl * pixelSize_.x()) < eps,
            photon.nscat() != 0 && std::abs(imagePosShifted.y() - yl * pixelSize_.y()) < eps,
            1.0);
    }
}


void Observer::bin(Photon const& photon, const Vector3d &pos1, const Vector3d &pos2)
{
    binMaskedLine(photon, pos1, pos2, 1.0);
}


void Observer::binHex(Photon const& photon, Vector3d const& pos1, Vector3d const& pos2)
{
    constexpr double a = 1.0 / 2.6457513110645907; // 1.0 / std::sqrt(7);
    constexpr double sqrt3 = 1.7320508075688772;

    constexpr Vector3d shifts[7] = {{0., 0., 0.}, {1.5*a, sqrt3/2.*a, 0}, {1.5*a, -sqrt3/2.*a, 0},
        {0, sqrt3*a, 0}, {0, -sqrt3*a, 0},
        {-1.5*a, sqrt3/2.*a, 0}, {-1.5*a, -sqrt3/2.*a, 0}};

    Vector3d const v = pos2 - pos1;
    double const r = v.norm() / 2.;
    Vector3d const x = Vector3d(-v.z(), v.x() * v.y() > 0 ? -v.z() : v.z(), v.x() * v.y() > 0 ? v.x() + v.y() : v.x() - v.y()).normalized();
    Vector3d const y = vectorProduct(v, x).normalized();

    for (int i=0; i!=7; ++i)
    {
        Vector3d const shift = r * (shifts[i].x() * x + shifts[i].y() * y);
        binMaskedLine(photon, pos1 + shift, pos2 + shift, 1.0 / 7.0);
    }
}


inline Vector2d Observer::project(Vector3d const &position) const
{
    return { position.y() * cosp_ - position.x() * sinp_,
        position.z() * direction_.sinTheta() - direction_.cosTheta() * (position.y()*sinp_ + position.x()*cosp_)};
}


inline void Observer::binToPixel(Photon const& photon, int64_t const x, int64_t const y, double const weight)
{
    result_.bin(photon, x, y, -1, weight);
    result0_.bin(photon, x, y, 0, weight);
    result1_.bin(photon, x, y, 1, weight);
    result2_.bin(photon, x, y, 2, weight);
}


inline void Observer::binPoint(const Photon &photon, int64_t const x, int64_t const y, bool const onXBorder, bool const onYBorder, double const weight)
{
    if (onXBorder && x > 0 && onYBorder && y > 0)
    {
        binToPixel(photon, x,   y,   0.25 * weight);
        binToPixel(photon, x,   y-1, 0.25 * weight);
        binToPixel(photon, x-1, y,   0.25 * weight);
        binToPixel(photon, x-1, y-1, 0.25 * weight);
    } else if (onXBorder && x > 0) {
        binToPixel(photon, x, y, 0.5 * weight);
        binToPixel(photon, x-1, y, 0.5 * weight);
    } else if (onYBorder && y > 0) {
        binToPixel(photon, x, y, 0.5 * weight);
        binToPixel(photon, x, y-1, 0.5 * weight);
    } else {
        binToPixel(photon, x, y, weight);
    }
}


inline void Observer::binLine(Photon const& photon, const Vector2d &pos1, const Vector2d &pos2, double const weight)
{
    double const x1 = nx_ * pos1.x() * rImageRev_;
    double const y1 = ny_ * pos1.y() * rImageRev_;
    auto const xl1 = static_cast<std::int64_t>(x1) - (x1 < 0.0);
    auto const yl1 = static_cast<std::int64_t>(y1) - (y1 < 0.0);

    double const x2 = nx_ * pos2.x() * rImageRev_;
    double const y2 = ny_ * pos2.y() * rImageRev_;
    auto const xl2 = static_cast<std::int64_t>(x2) - (x2 < 0.0);
    auto const yl2 = static_cast<std::int64_t>(y2) - (y2 < 0.0);

    double const eps = std::numeric_limits<float>::epsilon();

    if (xl1 == xl2 && yl1 == yl2)
    {
        binPoint(
            photon,
            xl2,
            yl2,
            std::abs(pos2.x() - pos1.x()) < eps && std::abs(pos2.x() - xl2 * pixelSize_.x()) < eps,
            std::abs(pos2.y() - pos1.y()) < eps && std::abs(pos2.y() - yl2 * pixelSize_.y()) < eps,
            weight);

        return;
    }

    double const dx = pos2.x() - pos1.x();
    double const dy = pos2.y() - pos1.y();
    int const xDir = dx > 0. ? 1.0 : -1.0;
    int const yDir = dy > 0. ? 1.0 : -1.0;

    int64_t borderX = xl1;
    int64_t const lastBorderX = xl2 + xDir;

    int64_t borderY = yl1;
    int64_t const lastBorderY = yl2 + yDir;

    double const totalW = std::sqrt(dx*dx + dy*dy);

    double x = pos1.x();
    double y = pos1.y();

    double xborder = dx > 0 ? std::min(static_cast<double>(borderX+1) * pixelSize_.x(), pos2.x())
            : std::max(static_cast<double>(borderX) * pixelSize_.x(), pos2.x());

    double yborder = dy > 0 ? std::min(static_cast<double>(borderY+1) * pixelSize_.y(), pos2.y())
            : std::max(static_cast<double>(borderY) * pixelSize_.y(), pos2.y());

    while (borderX != lastBorderX && borderY != lastBorderY)
    {
        double const xt = dx != 0 ? (xborder - x) / dx : std::numeric_limits<double>::infinity();
        double const yt = dy != 0 ? (yborder - y) / dy : std::numeric_limits<double>::infinity();

        if (xt < yt)
        {
            double yNew = y + xt * dy;
            double w = std::sqrt((yNew - y)*(yNew-y) + (xborder - x)*(xborder - x));

            binPoint(
                photon,
                borderX,
                borderY,
                false,
                dy < eps && std::fabs(y - borderY * pixelSize_.y()) < eps,
                weight * w / totalW);

            borderX += xDir;
            x = xborder;
            y = yNew;

            xborder = dx > 0 ? std::min(static_cast<double>(borderX+1) * pixelSize_.x(), pos2.x())
                    : std::max(static_cast<double>(borderX) * pixelSize_.x(), pos2.x());
        } else {
            double xNew = x + yt * dx;
            double w = std::sqrt((xNew - x)*(xNew-x) + (yborder - y)*(yborder - y));

            binPoint(
                photon,
                borderX,
                borderY,
                dx < eps && std::fabs(x - borderX * pixelSize_.x()) < eps,
                false,
                weight * w / totalW);

            borderY += yDir;
            x = xNew;
            y = yborder;

            yborder = dy > 0 ? std::min(static_cast<double>(borderY+1) * pixelSize_.y(), pos2.y())
                    : std::max(static_cast<double>(borderY) * pixelSize_.y(), pos2.y());
        }
    }
}


void Observer::binMaskedLine(Photon const& photon, const Vector3d &pos1, const Vector3d &pos2, double const lineWeight)
{
    Vector2d const imagePos1 = project(pos1);
    Vector2d const imagePos2 = project(pos2);
    Vector2d const shift{rImage_, rImage_};

    if (rMask2_ < std::numeric_limits<double>::epsilon())
    {
        binLine(photon, imagePos1 + shift, imagePos2 + shift, lineWeight);
        return;
    }

    if ((imagePos1 - imagePos2).norm2() < std::numeric_limits<double>::epsilon() && imagePos1.norm2() >= rMask2_)
    {
        binLine(photon, imagePos1 + shift, imagePos2 + shift, lineWeight);
        return;
    }

    double const r1sq = imagePos1.norm2();
    double const r2sq = imagePos2.norm2();

    if (r1sq <= rMask2_ && r2sq <= rMask2_)
    {
        return;
    }

    // line equation Ax + By + C = 0 based on two points
    double const A = imagePos1.y() - imagePos2.y();
    double const B = imagePos2.x() - imagePos1.x();
    double const C = imagePos1.x() * imagePos2.y() - imagePos2.x() * imagePos1.y();

    // closest point on the line
    double const A2B2 = A * A + B * B;
    Vector2d const closest{-A * C / A2B2, -B * C / A2B2};

    if (rMask2_ <= closest.norm2() + std::numeric_limits<double>::epsilon())
    {
        binLine(photon, imagePos1 + shift, imagePos2 + shift, lineWeight);
        return;
    }

    double const d2 = rMask2_ - C * C / A2B2;
    double const m = std::sqrt(d2 / A2B2);
    double const length = std::sqrt(A2B2);

    Vector2d const line{-B, A};

    double const cos1 = line * (closest - imagePos1);
    double const cos2 = line * (closest - imagePos2);

    if (cos1 * cos2 < 0)
    {
        double const cos1sign = (cos1 > 0.0) ? 1.0 : -1.0;
        // first subsegment
        if ((imagePos1 - closest).norm2() > d2)
        {
            Vector2d const border = closest + Vector2d{cos1sign * B * m, -cos1sign * A * m};
            binLine(photon, imagePos1 + shift, border + shift, lineWeight * (imagePos1 - border).norm() / length);
        }

        // second subsegment
        if ((imagePos2 - closest).norm2() > d2)
        {
            Vector2d const border = closest + Vector2d{-cos1sign * B * m, cos1sign * A * m};
            binLine(photon, imagePos2 + shift, border + shift, lineWeight * (imagePos2 - border).norm() / length);
        }
    }
    else
    {
        double const cosSign = (cos1 > 0.0 || cos2 > 0.0) ? 1.0 : -1.0;
        Vector2d const border = closest + Vector2d{cosSign * B * m, -cosSign * A * m};

        double const cos1b = line * (border - imagePos1);
        double const cos2b = line * (border - imagePos2);
        if (cos1b * cos2b < 0)
        {
            // intersects
            if (cos1b * cosSign > 0)
            {
                binLine(photon, imagePos1 + shift, border + shift, lineWeight * (imagePos1 - border).norm() / length);
            } else {
                binLine(photon, imagePos2 + shift, border + shift, lineWeight * (imagePos2 - border).norm() / length);
            }
        } else if (cos1b * cosSign >= 0 && cos2b * cosSign >= 0) {
            // outside
            binLine(photon, imagePos1 + shift, imagePos2 + shift, lineWeight);
        }
    }
}
