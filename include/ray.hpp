/// Global definitions
#define a 1.0       /// primary ray interaction probability scaling
#define F 0.5       /// primary ray interaction energy fraction
#define G 0.1       /// secondary ray deposition constant


class Ray
{
    public:
    /// Evolve ray by one step, return a tuple of:
    ///     - i,j for pixel visited
    ///     - energy deposited to that pixel
    ///     - new rays generated
    std::tuple<std::pair<int>, double, std::vector<Ray>> evolve();

    private:
    bool m_active;                          /// whether ray is active
    const bool m_primary;                   /// whether ray is primary
    const couble m_angle;                   /// angle of ray
    std::pair<int> m_current_pixel;         /// current pixel
    Simulation::PIXEL_EDGE m_current_edge;  /// current edge
    double m_current_edge_dist;             /// current distance along edge (from TOP or LEFT)

    /// Determine whether primary ray interacts at current pixel
    bool _interact();
    
    /// Trace ray through current pixel and return:
    ///     - distance travelled in pixel
    ///     - i,j for next pixel
    std::pair<double, std::pair<int>> _trace(); 

}
