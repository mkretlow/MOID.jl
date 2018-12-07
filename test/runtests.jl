using Test
using MOID
using Printf


@testset "Test MOID" begin

   SHOW = true

   rad = pi / 180

   elem = [ [" 1 Ceres"     , 2.7691652, 0.0760091,  73.59764,  80.30553,10.59407],
            ["29 Amphitrite", 2.5541136, 0.0726956,  63.36319, 356.34176, 6.08252],
            ["30 Urania"    , 2.3655722, 0.127581 ,  87.42605, 307.46872, 2.09575],
            ["50 Virginia"  , 2.6487939, 0.2859856, 200.08054, 173.52874, 2.83822],
            ["51 Nemausa"   , 2.3658354, 0.0675594,   2.58053, 175.9785 , 9.97718], ];

   # Expected results
   #Asteroid1  Asteroid2         MOID [AU]
   # 1 Ceres    1 Ceres          0.00000000000001
   # 1 Ceres   29 Amphitrite     0.15677463452737
   # 1 Ceres   30 Urania         0.24521440655832
   # 1 Ceres   50 Virginia       0.08934734026105
   # 1 Ceres   51 Nemausa        0.35972678460706

   result = [ 0.00000000000001, 0.15677463452737, 0.24521440655832, 0.08934734026105, 0.35972678460706 ]

   SHOW && println("\nAsteroid1  Asteroid2         MOID [AU]")

   for i in 1:length(elem)
      moid = moidF(elem[1][2:6]...,elem[i][2:6]...)
      @test moid ≈ result[i] atol=1E-14
      SHOW && @printf("%s   %-15s   %0.14f\n", elem[1][1], elem[i][1], moid)
   end

   SHOW && println()

end

# ---------------------------------------------------------------------------------------------------------------------------------#
