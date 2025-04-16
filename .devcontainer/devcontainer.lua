load("apptainer")

local img_name      = 'devcontainer_latest_amd64_2025_04_15_13_24_10.sif'
local img_directory = '/users/ataychameekiatchai/STAT-6358-Project/.devcontainer/' 
local img_path      = pathJoin(img_directory, img_name)

function build_command(cmd)
  local cmd_beginning = 'singularity exec '
  local cmd_ending    = img_path .. ' ' .. cmd
  local sh_ending     = ' "$@"'
  local csh_ending    = ' $*'
  local sh_cmd        = cmd_beginning .. cmd_ending .. sh_ending
  local csh_cmd       = cmd_beginning .. cmd_ending .. csh_ending
  set_shell_function(cmd, sh_cmd, csh_cmd)
end

-- build_command('bash')
