@testitem "Aqua" begin
    import Aqua
    Aqua.test_all(Wflow; persistent_tasks = false)
end
