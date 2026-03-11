@testitem "Aqua" begin
    import Aqua
    Aqua.test_all(WflowServer; persistent_tasks = false)
end
