
#ifndef CUSTOM_SCORE_TYPES_CUH
#define CUSTOM_SCORE_TYPES_CUH


struct ScoreType_u8x4{
    __host__ __device__
    operator unsigned int() const{
        return raw;
    }

    ScoreType_u8x4() = default;
    ScoreType_u8x4(const ScoreType_u8x4&) = default;

    __host__ __device__
    ScoreType_u8x4(unsigned int data){
        operator=(data);
    }

    __host__ __device__
    ScoreType_u8x4(std::uint8_t x, std::uint8_t y, std::uint8_t z, std::uint8_t w){
        set_x(x);
        set_y(y);
        set_z(z);
        set_w(w);
    }

    ScoreType_u8x4& operator=(const ScoreType_u8x4&) = default;

    __host__ __device__
    ScoreType_u8x4& operator=(unsigned int data){
        raw = data;
        return *this;
    }

    __host__ __device__
    std::uint8_t w() const{
        return (raw >> 24);
    }

    __host__ __device__
    std::uint8_t z() const{
        return (raw >> 16) & 0x000000FF;
    }

    __host__ __device__
    std::uint8_t y() const{
        return (raw >> 8) & 0x000000FF;
    }

    __host__ __device__
    std::uint8_t x() const{
        return (raw) & 0x000000FF;
    }

    __host__ __device__
    void set_w(std::uint8_t val){
        const unsigned int u = val;
        raw = (raw & 0x00FFFFFF) | (u << 24);
    }

    __host__ __device__
    void set_z(std::uint8_t val){
        const unsigned int u = val;
        raw = (raw & 0xFF00FFFF) | (u << 16);
    }

    __host__ __device__
    void set_y(std::uint8_t val){
        const unsigned int u = val;
        raw = (raw & 0xFFFF00FF) | (u << 8);
    }

    __host__ __device__
    void set_x(std::uint8_t val){
        const unsigned int u = val;
        raw = (raw & 0xFFFFFF00) | (u);
    }

    unsigned int raw{};
};

#endif