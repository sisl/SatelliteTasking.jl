push!(LOAD_PATH, "./")
push!(LOAD_PATH, "../")

# Julia Imports
using JSON
using Random

json_images = JSON.parsefile("./data/landsat_test.json")

num_images = length(keys(json_images["images"]))
num_samples = 3600

println("Found $num_images images")

sample_ids = Random.rand(1:16986, num_samples)

output_dict = Dict("images" => Any[])

for idx in Random.rand(1:16986, num_samples)
    push!(output_dict["images"], json_images["images"][idx])
end

JSON.print(open("./data/landsat_test_$num_samples.json", "w"), output_dict, 4)