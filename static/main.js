$(document).ready(function() {
            

    $('.search').autocomplete({
    source: function (request, response) {
        $.ajax({
            url: "/static/genes.json",
            dataType: 'json',
            data: request,
            success: function( data ) {
                 {
                    var filtered = data.filter(function (str) { 
                        return str.includes(request.term.toUpperCase())});
                    var aSearch = [];
                    // for each element in the main array ...
                    $(filtered).each(function(iIndex, sElement) {
                        // ... if element starts with input value ...
                        if (sElement.substr(0, request.term.length) == request.term.toUpperCase()) {
                            // ... add element
                            aSearch.push(sElement);
                        }
                    });

                    response(aSearch.splice(0, 6))
                };
            }
        }); 
       },  
       minLength: 1
      });   
});

