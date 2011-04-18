<?php
// This is an example of how to select shit from a database
/* > mysql> SELECT */
/* >     ->     product.product_info AS Product, */
/*    >     ->     product.product_price AS Price, */
/*    >     ->     SUM(orderrow.orderrow_quantity) AS Quantity, */
/*    >     ->     SUM(product.product_price*orderrow.orderrow_quantity) AS 'Total cost' */
/* >     -> FROM */
/*    >     ->     product, */
/*    >     ->     customer, */
/* >     ->     orderrow */
/* >     -> WHERE */
/* >     ->     product.product_id=orderrow.product_id */
/* >     -> AND */
/* >     ->     orderrow.customer_id=customer.customer_id */
/* >     -> AND */
/* >     ->     customer.customer_login='nobody' */
/* >     -> GROUP BY */
/*    >     ->     product.product_info; */
?>