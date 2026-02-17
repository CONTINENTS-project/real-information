#include <iostream>
#include <random>
#include <algorithm>
#include <array>
#include <bitset>
#include <iomanip>
#include "../bit_info.h"

#define n_elements 100
#define n_elem 100000
#define n_bits 32
int main() {
    float a[n_elem];
    float s[n_elem];
    double R[n_bits];

    auto print_a = [&a]() {
        for (int i = 0; i < 10; i++) std::cout << a[i] << ' ';
        std::cout << '\n';
    }; 
    auto print_s = [&s]() {
        for (int i = 0; i < 10; i++) std::cout << s[i] << ' ';
        std::cout << '\n';
    }; 
    auto print_R = [&R]() {
        for (auto elem : R) std::cout << elem << ' ';
        std::cout << '\n';
    }; 
    auto print_bits = [](std::bitset<n_bits> x) {
        for (int i = 0; i < n_bits; i++) std::cout << x[i] ;
        std::cout << '\n';
    };

    std::random_device rd;  // Will be used to obtain a seed for the random number engine
    std::mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()
    std::uniform_real_distribution<> dist(0.0, 1.0);
    for (int i = 0; i < n_elem; i++) {
        a[i] = dist(gen);
    }

    //std::sort(std::begin(a), std::end(a));
    //shave_template<float>(a, n_elem, 7, s);

    //redundancy(a, s, n_elem, R);
    //print_R();

    //double x = 12.148261, y;
    //shave_template<double>(&x, 1, 10, &y);
    //std::cout << "x " << x << " y " << y << '\n';

    //std::bitset<n_bits> mask = 0.15625;
    //print_bits(mask);

    //for (int i = 0; i < 32; i++) {
    //    shave_template<float>(a, n_elem, i, s);
    //    auto p = preserved_information_template<float>(a, s, n_elem);
    //    std::cout << "shaving " << i << '\n';
    //    print_a();
    //    print_s();
    //    std::cout << "preserved information " << p << '\n';
    //    std::cout << " " << '\n';
    //}

    //float x[n_elements] = { 0.18337965, 0.27502424, 0.6485108, 0.88652486, 0.05327642, 0.5636286, 0.9433723, 0.9619799, 
    //                        0.14971977, 0.5536389, 0.69815075, 0.76154107, 0.8933808, 0.76210606, 0.6540745, 0.25991088, 
    //                        0.17716247, 0.6412181, 0.076298475, 0.31624204};
    //float y[n_elements] = { 0.35556048, 0.7529169, 0.98594433, 0.9705406, 0.15364319, 0.57406586, 0.052820623, 0.7490818,
    //                        0.37220073, 0.63529676, 0.83195937, 0.80955184, 0.5724914, 0.5854158, 0.64927757, 0.62088704,
    //                        0.86874086, 0.15972632, 0.23554939, 0.6778932};

    float x[n_elements] = {1.07384448954147242E-4,  1.05271357656551255E-4,  1.04755144450627779E-4,  1.04055173171070791E-4,  1.03178324351598674E-4,  
                           1.02130912036995825E-4,  1.00927997572942483E-4,  9.9602221845461529E-5,  9.82359908815479257E-5,  9.69966512132148677E-5,  
                           9.62243431832830937E-5,  9.65566740526598407E-5,  9.75170612162669137E-5,  9.96829696619136826E-5,  -7.04904263596848192E-5,
                           3.54299356066262823E-4,  1.48545818340139611E-4,  5.94336928674805645E-5,  4.90843720509887582E-5,  5.71025994814323328E-5, 
                           5.93049570649650126E-5,  5.75760896470877258E-5,  5.35498194352349699E-5,  6.02067753733280631E-5,  6.47868059511549888E-5,  
                           3.11434422680530618E-5,  3.16632724318479327E-5,  3.79407667152461078E-5,  3.93113136711099162E-5,  3.78402310771351049E-5,  
                           3.5161806790606244E-5,  3.20057102556880945E-5,  2.85371472409434585E-5,  2.54491289590902121E-5,  2.36355545053613018E-5, 
                           2.08511986110449017E-5,  1.92826104620177545E-5,  2.03902792666211464E-5,  1.53157729707825669E-5,  1.15653177091072105E-5,  
                           9.08756068356392979E-6,  7.22431657569527101E-6,  5.66777237802826945E-6,  4.31690467059929054E-6,  3.13534359582819338E-6,  
                           2.12305282694963175E-6,  1.28458777988601084E-6,  6.27783466418706915E-7,  1.52998236678140362E-7,  -1.2549672730622003E-7,  
                           -2.09861874293941181E-7,  -8.27931884570093071E-8,  2.79621478258703529E-7,  9.42979067709946588E-7,  1.70693713380277403E-6,  
                           4.57516364356384377E-6,  1.90364486580216888E-6,  2.74495421430821084E-6,  7.29418668052367763E-6,  8.22008412742895018E-6,  
                           9.71206668106234871E-6,  1.18182772013005855E-5,  1.42001151093666046E-5,  1.67190074917682094E-5,  1.93315632916661174E-5,  
                           2.20178510897293822E-5,  2.47684167467732745E-5,  2.75645969350213841E-5,  3.03622906303843403E-5,  3.307688861504514E-5,  
                           3.55217180832312311E-5,  3.73803963221543663E-5,  3.73635723735900674E-5,  3.66570829453135322E-5,  1.69073031409335072E-5,  
                           6.33244700902557223E-5,  5.30997109915684774E-5,  -6.81987442292772991E-5,  3.56444939494101586E-6,  5.0365150876862473E-5,  
                           6.56253587228010121E-5,  7.16962827782111653E-5,  7.55726166400907289E-5,  7.88147476067385547E-5,  8.17507415710515088E-5,  
                           8.44749317468657427E-5,  8.70278305346413781E-5,  8.94278920768479073E-5,  9.16878475440019452E-5,  9.380163625252711E-5,  
                           9.57662105322596124E-5,  9.75809372073261207E-5,  9.92308765897652496E-5,  1.00708460229933426E-4,  1.02012879429872957E-4,  
                           1.03131127857672772E-4,  1.04061866037354371E-4,  1.0479682549490349E-4,  1.05340497107587046E-4,  1.07466891030908647E-4};
    float xs[n_elements];
    double data[32];
    shave_template(x, n_elements, 3, xs);
    bit_count_entropy(x, n_elements, data);
    redundancy(x, xs, n_elements, data);
    bitwise_real_information(x, n_elements, data);
    mutual_information(x, xs, n_elements, data);
    auto p = preserved_information_template<float>(x, xs, n_elements);

    float xr[n_elements];
    float xrs[n_elements];
    double data_r[32];
    //for (int i = 0; i < n_elements; i++) xr[i] = x[n_elements - i - 1];
    for (int i = 0; i < n_elements; i++) xr[i] = x[i];
    xr[0] = x[10];
    xr[10] = x[20];
    xr[20] = x[50];
    xr[50] = x[0];
    shave_template(xr, n_elements, 3, xrs);
    bit_count_entropy(xr, n_elements, data_r);
    redundancy(xr, xrs, n_elements, data_r);
    bitwise_real_information(xr, n_elements, data_r);
    mutual_information(xr, xrs, n_elements, data_r);
    auto pr = preserved_information_template<float>(xr, xrs, n_elements);

    std::cout << "mutual information\n" << data[0]  << " " << data_r[0]  << " " << data[0]  - data_r[0]  << "\n"
                                        << data[1]  << " " << data_r[1]  << " " << data[1]  - data_r[1]  << "\n"
                                        << data[2]  << " " << data_r[2]  << " " << data[2]  - data_r[2]  << "\n"
                                        << data[3]  << " " << data_r[3]  << " " << data[3]  - data_r[3]  << "\n"
                                        << data[4]  << " " << data_r[4]  << " " << data[4]  - data_r[4]  << "\n"
                                        << data[5]  << " " << data_r[5]  << " " << data[5]  - data_r[5]  << "\n"
                                        << data[6]  << " " << data_r[6]  << " " << data[6]  - data_r[6]  << "\n"
                                        << data[7]  << " " << data_r[7]  << " " << data[7]  - data_r[7]  << "\n"
                                        << data[8]  << " " << data_r[8]  << " " << data[8]  - data_r[8]  << "\n"
                                        << data[9]  << " " << data_r[9]  << " " << data[9]  - data_r[9]  << "\n"
                                        << data[10] << " " << data_r[10] << " " << data[10] - data_r[10] << "\n"
                                        << data[11] << " " << data_r[11] << " " << data[11] - data_r[11] << "\n"
                                        << data[12] << " " << data_r[12] << " " << data[12] - data_r[12] << "\n"
                                        << data[13] << " " << data_r[13] << " " << data[13] - data_r[13] << "\n"
                                        << data[14] << " " << data_r[14] << " " << data[14] - data_r[14] << "\n"
                                        << data[15] << " " << data_r[15] << " " << data[15] - data_r[15] << "\n"
                                        << data[16] << " " << data_r[16] << " " << data[16] - data_r[16] << "\n"
                                        << data[17] << " " << data_r[17] << " " << data[17] - data_r[17] << "\n"
                                        << data[18] << " " << data_r[18] << " " << data[18] - data_r[18] << "\n"
                                        << data[19] << " " << data_r[19] << " " << data[19] - data_r[19] << "\n"
                                        << data[20] << " " << data_r[20] << " " << data[20] - data_r[20] << "\n"
                                        << data[21] << " " << data_r[21] << " " << data[21] - data_r[21] << "\n"
                                        << data[22] << " " << data_r[22] << " " << data[22] - data_r[22] << "\n"
                                        << data[23] << " " << data_r[23] << " " << data[23] - data_r[23] << "\n"
                                        << data[24] << " " << data_r[24] << " " << data[24] - data_r[24] << "\n"
                                        << data[25] << " " << data_r[25] << " " << data[25] - data_r[25] << "\n"
                                        << data[26] << " " << data_r[26] << " " << data[26] - data_r[26] << "\n"
                                        << data[27] << " " << data_r[27] << " " << data[27] - data_r[27] << "\n"
                                        << data[28] << " " << data_r[28] << " " << data[28] - data_r[28] << "\n"
                                        << data[29] << " " << data_r[29] << " " << data[29] - data_r[29] << "\n"
                                        << data[30] << " " << data_r[30] << " " << data[30] - data_r[30] << "\n"
                                        << data[31] << " " << data_r[31] << " " << data[31] - data_r[31] << "\n";
    std::cout << "p " << p << '\n';
    std::cout << "pr " << pr << '\n';
    
    //for (int i = 0; i < 32; i++) {
    //    std::cout << data[i] << '\n';
    //}

    //for (auto elem : xs) std::cout << std::fixed << std::setprecision(8) << elem << '\n';

    pick_bits_to_shave_binary_search_template<float>(a, n_elem, 0.93, 0);

    return 0;
}