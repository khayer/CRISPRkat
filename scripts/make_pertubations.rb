require "csv"

# ruby scripts/make_pertubations.rb treatment treatment2 \
# 04quant/mageck.count.txt treatment_rep_1,treatment_rep_2 \
# treatment2_rep_1,treatment2_rep_2 > 05norm/treatment_treatment2_pertubation.txt

cond1 = ARGV[3].split(",")
cond2 = ARGV[4].split(",")

puts "Barcode Sequence\tScore"
kat = {}

CSV.foreach(ARGV[2], { :col_sep => "\t" , :headers => :first_row }) do |row|
	#STDERR.puts "KAT"
	cond1_elements = []
	cond2_elements = []
	cond1.each {|e| cond1_elements << row[e]}
	cond2.each {|e| cond2_elements << row[e]}
  average_cond1 = cond1_elements.inject{ |sum, el| sum + el }.to_f / cond1_elements.size
  #STDERR.puts average_cond1
	average_cond2 = cond2_elements.inject{ |sum, el| sum + el }.to_f / cond2_elements.size
	#STDERR.puts average_cond2
	fold_change = (average_cond1+1.0)/(average_cond2+1.0)
	#STDERR.puts fold_change
	kat[row["sgRNA"]]  = fold_change
	#puts "#{row["sgRNA"]}\t#{fold_change}"
end

#STDERR.puts kat
kat = kat.sort_by {|k,v| v}

### SORTING IS MISSING
kat.each do |k|
	puts "#{k[0]}\t#{k[1]}"
end