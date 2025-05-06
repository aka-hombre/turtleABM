from mesa.batchrunner import batch_run
import pandas as pd

from models.movingmodel import MovingModel

if __name__ == "__main__":
    # Testing Cases
    cases = [
        (0,   360),
        (0,   180),
        (180, 360),
        (45,  135),
        (225, 315),
        (330, 360),
    ]
    
    # Defining seeds for reproducibility
    seeds = list(range(42, 52)) 

    all_results = []

    # Loop through each deg_min, deg_max pair
    for deg_min, deg_max in cases:
        params = {
            "num_agents": [1],    # One agent per run
            "deg_min": [deg_min],  
            "deg_max": [deg_max],  
            "seed": seeds,        
        }

        # Run each configuration with 1 iteration per seed
        results = batch_run(
            model_cls=MovingModel,
            parameters=params,
            iterations=3,              
            max_steps=150*24,             
            number_processes=None,     
            data_collection_period=1,   
            display_progress=True,     # Show progress
        )
        
        all_results.extend(results)  # Collect results

    df = pd.DataFrame(all_results)
    df.to_csv("expiriment_results.csv", index=False)  
    print("Done: results in deg_cases_results.csv")
