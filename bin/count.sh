query="$1"
filename="$2"

function usage() {
  echo "usage: $0 <query> <filename>

counts the number of instances of <query> in <filename> via grep"

  exit 1
}

if [ "$query" == "" ]; then
  usage
fi

if [ "$filename" == "" ]; then
  usage
fi

echo "$query" "$(basename $filename)" "$(grep "$query" "$filename" | wc -l)"
