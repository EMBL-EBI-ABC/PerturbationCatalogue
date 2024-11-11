import { Component, OnInit } from '@angular/core';
import { HeaderComponent } from "../shared/header/header.component";
import { MatCardModule } from "@angular/material/card";
import { FooterComponent } from "../shared/footer/footer.component";
import { ElasticService } from '../elastic.service';
import { MatTableModule } from '@angular/material/table';
import { MatPaginator, MatPaginatorModule, PageEvent } from '@angular/material/paginator';
import { CommonModule } from '@angular/common';




@Component({
  selector: 'app-home',
  standalone: true,
  imports: [HeaderComponent, FooterComponent, MatCardModule, MatTableModule, MatPaginator, MatPaginatorModule, CommonModule,
],
  templateUrl: './home.component.html',
  styleUrl: './home.component.css'
})
export class HomeComponent implements OnInit {
  data: any;
  totalResults = 0;
  pageSize = 10;
  currentPage = 0;

  constructor(private elasticService: ElasticService) {}

  ngOnInit() {
    this.fetchData();
  }

  getKeys(obj: any): string[] {
    return Object.keys(obj);
  }

  fetchData() {
    const start = this.currentPage * this.pageSize;
    const size = this.pageSize;
    this.elasticService.getData(start, size).subscribe(
      (response) => {
        this.data = response.results;
        this.totalResults = response.total;
      },
      (error) => {
        console.error('Error fetching data:', error);
      }
    );
  }

  onPageChange(event: PageEvent) {
    this.pageSize = event.pageSize;
    this.currentPage = event.pageIndex;
    this.fetchData();
  }

}
