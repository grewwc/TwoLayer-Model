const exec = require('child_process').exec;

const all_stars = ['j0218', 'j1939', 'b1821'];

(function auto() {
    if (process.argv.length !== 3) {
        throw "you should have a parameter";
    }
    let star_name = process.argv[2];
    if (all_stars.includes(star_name)) {
        exec(`./main model_${star_name}`, (err, stdout, stderr) => {
            console.log(stdout, '\n', stderr);
            exec(`python /home/wwc129/work/python/${star_name}.py`, (err, stdout, stderr) => {
                console.log(stdout, '\n', stderr);
            })
        });
    } else {
        let dash = star_name.indexOf("_");
        let first_part = star_name.slice(dash+ 1);
        let second_part = star_name.slice(0, dash);
        let outer_operation;
        let changed_starName = star_name.slice(dash + 1) + "_" + star_name.slice(0, dash);
        if(second_part === "cur"){
            outer_operation = "cur";
        }else{
            outer_operation = changed_starName;
        }
        
        exec(`./main ${outer_operation}`, (err, stdout, stderr) => {
            console.log(stdout, '\n', stderr);
            if (dash === -1) {
                throw "naming convention is wrong";
            } else {
                exec(`/home/wwc129/Twolayer/data/${changed_starName}.sh`, (err, stdout, stderr) => {
                    console.log(stdout, '\n', stderr);
                })
            }
        });
    }
})();